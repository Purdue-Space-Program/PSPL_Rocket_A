from math import ceil, sqrt, atan, pi, degrees, floor
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import vehicle_parameters as vehicle


LB2KG = 0.453592
FT2M = 0.3048
IN2M = 0.0254
LBF2N = 4.44822
MPH2MPS = 0.44704

# Need q = 1/2 rho v^2 (done)
# Need AOA = tan^-1(wind_gust / velocity) (done)
# Need S (cross-sectional area)= pi * (diameter / 2)^2 (done)
# Need SD stability derivative (done)
# SD for fins = Cambridge equation 31 (done)
# Need Kfb for finSD = Cambridge equation 32 (done)
# SD for nose = 2 * machcoeff (done)
# machcoeff = 1 / sqrt(abs(mach^2 - 1)) (done)
# Need lift = q * S * AOA * SD
# Need total mass
# Need total length
# Need center of gravity = (sum mass_x * x_coord for all x in length) / total mass
# Need I rotational inertia = sum mass_x * (x_coord - center of gravity)^2
# Need ay lateral acceleration = sum lift / total mass
# Need r angular acceleration = (sum lift_x * (cp_x_coord - center of gravity)) / rotational inertia
# Need cp center of pressure
# Need finCP = Cambridge equation 33
# Need noseCP = noseLength / 2 literally middle of nose
# Need shear = (sum lift) - ay * mass up to x - r * mass up to x * (cg - mass up to x) for all x in length
# Need axial = - thrust + nose drag + fuselage drag + base drag + S * pressure difference + forward acceleration * mass up to x
# Given thrust, forward acceleration
# Worry about axial later

# Need rho, velocity, wind gust, diameter, fin dimensions, mach, linear density, total mass, total length, 

# Calculate dynamic pressure
def calcQ(rho, velocity, wind_gust):
    '''
    rho: air density [kg / m^3]
    velocity: air velocity = rocket velocity [m/s]
    Q: Dynamic pressure [Pa]
    '''
    return rho * (sqrt(velocity**2 + wind_gust**2))**2 / 2

# Calculate angle of attack
def calcAOA(wind_gust, velocity):
    '''
    wind_gust: wind gust velocity [m/s]
    velocity: rocket velocity [m/s]
    AOA: Angle of attack [radians]
    '''
    return atan(wind_gust / velocity)

# Calculate cross sectional area
def calcS(diameter):
    '''
    diameter: rocket diameter [m]
    S: Cross sectional area [m^2]
    '''
    return pi * (diameter / 2)**2

# Calculate fin stability derivative
def calcFinSD(root_chord, tip_chord, sweep_length, fin_height, numFins, diameter):
    '''
    root_chord: fin root chord length [m]
    tip_chord: fin tip chord length [m]
    sweep_length: fin sweep length [m]
    fin_height: fin height aka fin span [m]
    numFins: Number of fins [unitless]
    diameter: Rocket diameter [m]
    finSD: Stability derivative [unitless]
    '''
    mid_chord = sqrt(fin_height**2 + ((tip_chord - root_chord) / 2 + sweep_length)**2) # [m] Length of fin mid-chord line, Rocket Fin Design equation 1
    print(f"Fin mid-chord: {mid_chord}") # TEST
    Kfb = 1 + (diameter / 2) / (fin_height + diameter / 2) # Coefficient, Cambridge equation 32
    print(f"Fin Kfb: {Kfb}") # TEST
    stability_derivative = Kfb * (4 * numFins * (fin_height / diameter)**2) / (1 + sqrt(1 + (2 * mid_chord / (root_chord + tip_chord))**2)) # Fin stability derivative, Cambrdige Aerodynamic Equations equation 31
    return stability_derivative

# Calculate mach coefficient
def calcMachCoeff(coeff, mach):
    '''
    coeff: [unitless] Some coefficient, I think we use 1
    mach: [unitless] Speed of rocket in mach, something about the Prandtl-Glauert rule
    '''
    return coeff / sqrt(1 - mach**2) # Cambridge equation 55, 56, assume rocket goes subsonic

# Calculate nose stability derivative
def calcNoseSD(machCoeff):
    '''
    machCoeff: Mach coefficient
    why don't I include the mach coefficient for fins?
    '''
    return 2 * machCoeff # Cambridge equation 25

# Calculate lift
def calcLift(Q, S, AOA, SD):
    '''
    Q: Dynamic pressure [Pa]
    S: Cross sectional area [m^2]
    AOA: Angle of attack [radians]
    SD: Stability derivative
    Lift: Lift [N]
    '''
    return Q * S * AOA * SD # [N] Aspire page 14

# Calculate center of gravity
def calcCG(linear_density_array, length_along_rocket_linspace):
    '''
    linear_density_array: Array of linear density as a function of length [kg / m]
    length_along_rocket_linspace: Array of length along rocket [m]
    cg: Location of center of gravity of rocket from aft [m]
    '''
    dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]
    totalMass = np.sum(linear_density_array * dx)
    # print(totalMass / LB2KG)
    lengths = np.array(length_along_rocket_linspace)
    masses = np.array(linear_density_array * dx)
    moments = np.sum(lengths * masses)
    cg = moments / totalMass
    return cg

# print(calcCG(vehicle.linear_density_array, vehicle.length_along_rocket_linspace) / FT2M) # TEST

# Create an array of center of gravity with respect to time
def updateCG(vehicle, burn_time, total_time):
    '''
    vehicle: imported vehicle
    burn_time: burn time [s]
    '''
    length_along_rocket_linspace = vehicle.length_along_rocket_linspace
    cg = []
    dt = 0.005
    times = np.arange(0.0, total_time, dt)
    mass_flow_rate = vehicle.parameters.mass_flow_rate
    OF_ratio = vehicle.parameters.OF_ratio
    ox_flow_rate = mass_flow_rate * OF_ratio / (1 + OF_ratio)
    fuel_flow_rate = mass_flow_rate * 1 / (1 + OF_ratio)
    ox_location, ox_length = vehicle.oxidizer_tank.bottom_distance_from_aft, vehicle.oxidizer_tank.length
    fuel_location, fuel_length = vehicle.fuel_tank.bottom_distance_from_aft, vehicle.fuel_tank.length
    # print(times) # TEST
    for t in times:
        if t > burn_time: 
            cg.append(cg[-1])
        else:
            linear_density_array = np.zeros(vehicle.num_points)
            for component in vehicle.mass_distribution.components:
                mass = component.mass
                if component.name == 'oxidizer_tank':
                    mass -= ox_flow_rate * t
                if component.name == 'fuel_tank':
                    mass -= fuel_flow_rate * t
                linear_density = mass / component.length
                for index, length_along_rocket in enumerate(length_along_rocket_linspace):
                    above_component_bottom = length_along_rocket >= component.bottom_distance_from_aft
                    below_component_top = length_along_rocket <= (component.bottom_distance_from_aft + component.length)
                    
                    if (above_component_bottom and below_component_top):
                        linear_density_array[index] += linear_density
            cg.append(calcCG(linear_density_array, length_along_rocket_linspace))
    return cg

'''
TEST
print(updateCG(vehicle, 2.24, 10)[0]) # TEST & VALIDATE
print(calcCG(vehicle.linear_density_array, vehicle.length_along_rocket_linspace))
plt.plot(np.arange(0.0, 10, 0.005), updateCG(vehicle, 2.24, 10))
plt.xlabel("time")
plt.ylabel("center of gravity")
plt.show()
'''

# Calculate rotational inertia
def calcRotationalInertia(linear_density_array, length_along_rocket_linspace, cg):
    '''
    linear_density_array: Array of linear density across rocket [array]
    length_along_rocket_linspace: Numpy linspace for lengths along rocket [array]
    cg: Location of center of gravity [m]
    '''
    dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]
    mass_model = linear_density_array * dx
    inertia = 0
    for x in range(len(length_along_rocket_linspace)):
        inertia += mass_model[x] * (length_along_rocket_linspace[x] - cg)**2
    return inertia

# Calculate lateral acceleration, need a way to update total_mass
def calcLateralAcceleration(noseLift, finLift, total_mass):
    '''
    lift_dict: Dictionary of lift forces
    total_mass: Total mass of rocket
    '''
    return (noseLift + finLift) / total_mass # [m/s^2] a = F / m

# Calculate the location of the fin center of pressure
def calcFinCP(root_chord, tip_chord, sweep_length, fin_height, total_length, noseconeToFin):
    '''
    root_chord: Length of fin root chord [m]
    tip_chord: Length of fin tip_chord [m]
    sweep_length: Length of fin sweep length [m]
    fin_height: Fin height [m]
    total_length: Total length of rocket [m]
    noseconeToFin: Length from noescone to find [m]
    finCP: Location of center of pressure of the fin
    '''
    mid_chord = sqrt(fin_height**2 + ((tip_chord - root_chord) / 2 + sweep_length)**2) # [m] Length of fin mid-chord line, Rocket Fin Design equation 1
    first_term = noseconeToFin
    second_term = (mid_chord / 3) * ((root_chord + 2 * tip_chord) / (root_chord + tip_chord))
    third_term = (1 / 6) * (root_chord + tip_chord - (root_chord * tip_chord) / (root_chord + tip_chord))
    
    finCP = first_term + second_term + third_term # Cambridge equation 33
    # print(length) # TEST
    # print(cp) # TEST
    return total_length - finCP # [m] Fin CP from aft

# Calculate the location of the nose center of pressure
def calcNoseCP(nosecone_length, total_length):
    '''
    nosecone_length: Length of nosecone [m]
    total_length [m]
    noseCP: Location of center of pressure of the nose
    '''
    noseCP = 0.5 * nosecone_length # [m] Nosecone CP, Cambridge equation 28
    return total_length - noseCP # [m] Nose CP from aft

# Calculate angular acceleration across rocket length
def calcAngularAcceleration(noseLift, finLift, noseCP, finCP, inertia, cg):
    '''
    noseLift: Nosecone lift [N]
    finLift: Fin lift [N]
    noseCP: Location of center of pressure of the nose [m]
    finCP: Location of center of pressure of the fin [m]
    inertia: Rotational inertia around center of gravity [kg m^2]
    cg: Location of center of gravity of rocket [m]
    r: Angular acceleration [1 / s^2]
    '''
    r = ((-1) * noseLift * (abs(noseCP - cg)) + finLift * (abs(finCP - cg))) / inertia
    return r

# Calculate shear forces across the rocket length and output an array of shear forces
def calcShear(noseLift, finLift, noseCP, finCP, ay, linear_density_array, length_along_rocket_linspace, r, cg):
    '''
    noseLift: Nosecone lift [N]
    finLift: Fin lift [N]
    noseCP: Location of center of pressure of the nose [m]
    finCP: Location of center of pressure of the fin [m]
    ay: Lateral acceleration [m / s^2]
    linear_density_array: Array of linear density across rocket length [kg / m]
    length_along_rocket_linspace: Array of rocket lengths [m]
    shear_array: Array of shear forces across rocket length [N]
    '''
    dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]
    cg_rel_lengths = np.array(-1 * length_along_rocket_linspace + cg)
    mass_model = np.cumsum(linear_density_array * dx) # aft to nose
    cumulative_moment_about_cg = np.cumsum(linear_density_array * dx * cg_rel_lengths) # aft to nose
    shear_array = (-1) * ay * mass_model - r * cumulative_moment_about_cg
    #print(np.cumsum(linear_density_array * dx))
    # print(shear_array) # TEST
    #print(f"Lateral acceleration: {ay}")
    #print(f"Nose lift: {noseLift}") # TEST
    #print(f"Fin lift: {finLift}") # TEST
    shear_array[int(noseCP / dx):] += noseLift
    shear_array[int(finCP / dx):] += finLift

    return shear_array

# Calculate bending forces across rocket length and output an array of bending forces
def calcBending(shear_array, length_along_rocket_linspace):
    '''
    shear_array: Array of shear forces across rocket length [N]
    length_along_rocket_linspace: Array of rocket lengths [m]
    bending_array: Array of bending forces across rocket length [N m]
    '''
    dy = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]
    bending_array = np.cumsum(shear_array) * dy
    return bending_array

# Calculate axial forces across rocket length and output an array of bending forces
def calcAxial(thrust, ax, linear_density_array, length_along_rocket_linspace, rho, cd, S, velocity):
    '''
    thrust: Rocket thrust [N]
    ax: Acceleration [m / s^2]
    linear_density_array: Array of linear density across rocket length [kg / m]
    length_along_rocket_linspace: Array of rocket length [m]
    axial: Array of axial forces across rocket length [N]
    '''
    dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]
    masses = np.array(linear_density_array) * dx # aft to nose
    masses1 = masses[::-1] # nose to aft
    masses2 = np.cumsum(masses1) # nose to aft
    masses3 = masses2[::-1] # aft to nose

    axial_tb = []
    for mass in masses3:
        axial_p = (ax * mass) + mass * 9.81
        axial_tb.append(axial_p)
    fin_drag = 0.5 * rho * cd * S * velocity**2 # Fin drag [N]
    axial_tb = np.array(axial_tb) + fin_drag # aft to nose
    axial = axial_tb # aft to nose
    # to_strut = floor(32.2 * IN2M / dx)
    # print(to_strut)
    # print(axial[to_strut])
    return axial

'''
plt.plot(vehicle.length_along_rocket_linspace, calcAxial(vehicle.parameters.jet_thrust, vehicle.parameters.max_acceleration, vehicle.linear_density_array, vehicle.length_along_rocket_linspace))
plt.xlabel("length [m]")
plt.ylabel("axial force [N]")
plt.show()
'''

'''
def codeFriendlyName(str, delimiter='_'):
    return delimiter.join(str.split(' (')[0].split(' ')).lower() # Delimiter is character between words ex. "_"

def getRocketDict(vehicle):
    chamber_pressure = vehicle.parameters.chamber_pressure # [Pa]
    thrust = vehicle.parameters.jet_thrust # [N]
    tank_pressure = vehicle.parameters.tank_pressure # [Pa]
    fuel_tank_length = vehicle.parameters.fuel_tank_length # [m]
    fuel_total_mass = vehicle.parameters.fuel_total_mass # [kg]
    oxidizer_tank_length = vehicle.parameters.oxidizer_tank_length # [m]
    oxidizer_total_mass = vehicle.parameters.oxidizer_total_mass # [kg]
    total_length = vehicle.parameters.total_length # [m]
    wet_mass = vehicle.parameters.wet_mass # [kg]
    dry_mass = vehicle.parameters.dry_mass # [kg]
    off_the_rail_acceleration = vehicle.parameters.off_the_rail_acceleration * 9.81 # [m/s^2] obtained by multiplying by gravitational acceleration
    off_the_rail_velocity = vehicle.parameters.off_the_rail_velocity # [m/s]
    max_acceleration = vehicle.parameters.max_acceleration * 9.81 # [m/s^2] obtained by multiplying by gravitational acceleration
    max_velocity = vehicle.parameters.max_velocity # [m/s]

    lox_tank = {'mass': 0, 'length': oxidizer_tank_length, 'prop_mass': oxidizer_total_mass}
    fuel_tank = {'mass': 0, 'length': fuel_tank_length, 'prop_mass': fuel_total_mass}
    

# Unpack 'Vehicle Sections' in rocket Excel file into a dictionary
def getRocketSections(xlsx):
    rocket_sections = {}
    rocket_inputs = pd.read_excel(xlsx, sheet_name='Vehicle Sections') # Opens and reads a sheet called "Vehicle Sections"

    for section in rocket_inputs.iloc: # Iterates through each row (rocket section) of the sheet
        sec_name = codeFriendlyName(section['Name'])
        # if sec_name == 'nosecone': continue # Nosecone is ejected in recovery
        # print(sec_name) # TEST
        section_dict = section.drop('Name').to_dict() # Each rocket section will have the format { mass: xx, length: xx, ... }

        rocket_sections[sec_name] = {codeFriendlyName(k):v for k,v in section_dict.items() if v==v} # Creates dictionary for each rocket section, eliminates dictionary items wihout values (empty cells in spreadsheet)
        # Convert units to metric
        if 'mass' in rocket_sections[sec_name]:
            rocket_sections[sec_name]['mass'] *= LB2KG
        if 'prop_mass' in rocket_sections[sec_name]:
            rocket_sections[sec_name]['mass'] += 0.08 * rocket_sections[sec_name]['prop_mass'] * LB2KG # Need to know what to change 0.08 to and how much prop is left at max q
        if 'length' in rocket_sections[sec_name]:
            rocket_sections[sec_name]['length'] *= FT2M
        if 'thickness' in rocket_sections[sec_name]:
            rocket_sections[sec_name]['thickness'] *= IN2M

    return rocket_sections # [dictionary]

# Unpack 'Point Masses' in rocket Excel file into a dictionary
def getPointMasses(xlsx):
    point_mass_inputs = pd.read_excel(xlsx, sheet_name='Point Masses', header=0, usecols=[1,2])
    point_masses = []
    for point_mass in point_mass_inputs.iloc: # Iterates through each row (point mass) of the sheet
        mass, xcoord = list(point_mass)
        point_masses.append({'mass':mass * LB2KG, 'x_coordinate':xcoord * FT2M}) # Only point mass is batteries
    # print(point_masses) # TEST
    return point_masses # [dictionary]

# Unpack 'Aerodynamic Properties' in rocket Excel file into a dictionary
def getAeroProperties(xlsx):
    location_properties = {}
    location_inputs = pd.read_excel(xlsx, sheet_name='Aerodynamic Properties') # Opens and reads a sheet called "Aerodynamic Properties"

    for location in location_inputs.iloc: # Iterates through each row (location) of the sheet
        location_name = codeFriendlyName(location['Name'])
        # print(location_name) # TEST
        location_dict = location.drop('Name').to_dict() # Each location will have the format { velocity: xx, acceleration: xx, ... }

        location_properties[location_name] = {codeFriendlyName(k):v for k,v in location_dict.items()} # Creates a dictionary for each location
        # Convert units to metric
        if 'thrust' in location_properties[location_name]:
            location_properties[location_name]['thrust'] *= LBF2N
        if 'max_wind_guss' in location_properties[location_name]:
            location_properties[location_name]['max_wind_gust'] *= MPH2MPS
    # print(location_properties) # TEST
    return location_properties # [dictionary]

# Create a mass model of rocket in the form of an array where each point in the array has a mass of the rocket at that location
def getMassModel(rocket_dict, point_masses, dy):
    element_length = dy # Set element length of 5 mm
    mass_model = []
    model_index = 0
    pointMassIndices = [ceil(pm['x_coordinate'] / element_length) for pm in point_masses] # Find location in mass array that corresponds to point mass x-coordinate on rocket
    # print(pointMassIndices) # TEST
    for section in reversed(rocket_dict.values()): # Go fin to nose instead of nose to fin
        secMass = section.get('mass')
        secLength = section.get('length')
        if not secMass: continue # Skips fins
        elementsPerSection = ceil(secLength / element_length)
        massPerElement = secMass / elementsPerSection # Split a section mass across the entire section
        for i in range(elementsPerSection):
            if model_index in pointMassIndices: # Add mass of point masses to mass_model
                point_mass_mass = point_masses[pointMassIndices.index(model_index)]['mass']
                mass_model.append(massPerElement + point_mass_mass)
            else:
                mass_model.append(massPerElement) # Add the element mass to mass_model
            model_index += 1
    # print(mass_model) # TEST
    return np.array(mass_model) # [array]

# Calculate total mass of rocket
def getTotalMass(rocket_dict, point_masses):
    totalPointMass = sum([pm['mass'] for pm in point_masses]) # Sum point masses
    totalSectionMass = sum([rocket_section.get('mass', 0) for rocket_section in rocket_dict.values()]) # Sum section masses
    return totalPointMass + totalSectionMass # [kg] Sum point masses and section masses

# Calculate total length of rocket
def getTotalLength(rocket_dict):
    totalLength = sum([rocket_section.get('length', 0) for rocket_section in rocket_dict.values()]) # Sum section lengths
    return totalLength # [m]

# Calculate center of gravity of rocket
def getCG(rocket_dict, point_masses, totalMass):
    moment = 0
    length = 0
    for section in reversed(rocket_dict.values()): # Go from fin to nose instead of nose to fin
        secMass = section.get('mass', 0)
        secLength = section.get('mass', 0)
        xCoord = length + secLength / 2 # Each section can be represented as a point mass
        moment += secMass * xCoord # Multiply masses by x coordinate to get moment for each section
        length += secLength # Go to next section
    
    pointMass = np.array([pm['mass'] for pm in point_masses])
    pointCoord = np.array([pm['x_coordinate'] for pm in point_masses])
    # Multiply masses by x coordinate to get moment for each point mass
    pointMoments = pointMass * pointCoord
    moment += np.sum(pointMoments)

    cg = moment / totalMass # Center of gravity calculation
    return cg # [m]

# Fin stability derivative for trapezoidal fins
def getFinSD(rocket_dict, rocket_diameter):
    root_chord, tip_chord, sweep_length, fin_height, n = rocket_dict['fins']['root_chord'], rocket_dict['fins']['tip_chord'], rocket_dict['fins']['sweep_length'], rocket_dict['fins']['height'], rocket_dict['fins']['number_of_fins'] # Unpack rocket_dict
    mid_chord = sqrt(fin_height**2 + ((tip_chord - root_chord) / 2 + sweep_length)**2) # [m] Length of fin mid-chord line, Rocket Fin Design equation 1
    Kfb = 1 + (rocket_diameter / 2) / (fin_height + rocket_diameter / 2) # Coefficient, Cambridge equation 32

    stability_derivative = Kfb * (4 * n * (fin_height / rocket_diameter)**2) / (1 + sqrt(1 + (2 * mid_chord / (root_chord + tip_chord))**2)) # Fin stability derivative, Cambrdige Aerodynamic Equations equation 31
    return stability_derivative

# Fin center of pressure from Barrowman Equation
def getFinCP(noseconeToFin, rocket_dict, length):
    root_chord, tip_chord, sweep_length, fin_height = rocket_dict['fins']['root_chord'], rocket_dict['fins']['tip_chord'], rocket_dict['fins']['sweep_length'], rocket_dict['fins']['height'] # Unpack rocket_dict
    mid_chord = sqrt(fin_height**2 + ((tip_chord - root_chord) / 2 + sweep_length)**2) # [m] Length of fin mid-chord line, Rocket Fin Design equation 1
    first_term = noseconeToFin
    second_term = (mid_chord / 3) * ((root_chord + 2 * tip_chord) / (root_chord + tip_chord))
    third_term = (1 / 6) * (root_chord + tip_chord - (root_chord * tip_chord) / (root_chord + tip_chord))

    cp = first_term + second_term + third_term # Cambridge equation 33
    # print(length) # TEST
    # print(cp) # TEST
    return length - cp # [m] Fin CP from aft

# Calculate corrected aerodynamic coefficient
def getMachAdjustedCoeff(coeff, mach):
    return coeff / sqrt(abs(mach**2 - 1)) # Cambridge equation 55, 56

# Calculate compressible flow stability derivative for conical, ogive, and parabolic nosecone
def getNoseSD(machCoeff):
    return 2 * machCoeff # Cambridge equation 25

# Calculate center of pressure of nosecone
def getNoseCP(nosecone_length, length): # Parabolic nosecone CP equation
    xCoord = length - 0.5 * nosecone_length # [m] Nosecone CP from aft, Cambridge equation 28
    return xCoord

# Calculate the total mass of the rocket
def getTotalMass(rocket_dict):
    totalMass = 0
    for sec_name in rocket_dict: # Loop through rocket_dict and sum masses
        if 'mass' in rocket_dict[sec_name]:
            totalMass += rocket_dict[sec_name]['mass'] # [kg]
    return totalMass

# Calculate angle of attack when rocket encounters a wind gust
def getAOA(aero_dict, location):
    return degrees(atan(aero_dict[location]['max_wind_gust'] / aero_dict[location]['velocity'])) # [degrees]

# Calculate dynamic pressure
def getQ(aero_dict, location):
    velocity, rho = aero_dict[location]['velocity'], aero_dict[location]['air_density']
    return rho * velocity**2 / 2 # [N / m^2]

# Calculate cross-sectional area
def getArea(diameter):
    return pi * (diameter / 2)**2 # [m^2]

# Calculate lift force
def getLiftForce(Q, S, AOA, SD):
    return Q * S * AOA * SD # [N] Aspire page 14

# Calculate rotational inertia
def getRotationalInertia(mass_model, cg, length):
    lengths = np.linspace(0, length, len(mass_model))
    rotationalInertia = np.sum(mass_model * (lengths - cg)**2)
    return rotationalInertia # [kg m^2]

# Calculate lateral acceleration when rocket encounters wind gust
def getLatAccel(lift_dict, totalMass):
    noseLift, finLift, boattailLift = lift_dict['nose'], lift_dict['fin'], lift_dict['boattail'] # Unpack lift_dict
    return (noseLift + finLift + boattailLift) / totalMass # [m/s^2] a = F / m

# Calculate angular acceleration when rocket encounters wind gust
def getAngularAccel(lift_dict, cp_dict, cg, inertia):
    lift = np.array([v for k,v in lift_dict.items()]) # Create an array for lift from lift_dict
    cp = np.array([v for k,v in cp_dict.items()]) # Create an array for cp from cp_dict
    cpLocation = np.absolute(cp - cg) # Create an array of distance of cp from cg
    sum = 0
    for i in range(len(lift_dict) - 1):
        sum += (-1)**(i + 1) * lift[i] * cpLocation[i] # lift_dict should be nose, fin, boattail so nose (-), fin (+), boattail (-)
    return sum / inertia # [1 / s^2]

# Calculate shear force
def getShearForce(mass_model, totalLength, lift_dict, cp_dict, ay, r, cg):
    dy = totalLength / len(mass_model)
    lengths = np.linspace(0, totalLength, len(mass_model))
    # print(lengths) # TEST
    
    # Calculate shear force at each point on rocket Aspire page 21
    shear_array = (-1) * (ay * np.cumsum(mass_model))# + r * np.cumsum(mass_model * (cg - lengths)))
    shear_array[int(cp_dict['nose'] / dy):] += lift_dict['nose']
    shear_array[int(cp_dict['fin'] / dy):] += lift_dict['fin']
    shear_array[int(cp_dict['boattail'] / dy):] -= lift_dict['boattail']
    
    # print(shear_array) # TEST
    return shear_array # [array, N]

# Calculate bending force
def getBendingForce(shear_array, totalLength):
    dy = totalLength / len(shear_array)
    bending_array = np.cumsum(shear_array) * dy # Bending graph is integral of shear graph
    return bending_array # [array, N * m]

def getAxialForces(aero_dict, mass_model, S):
    velocity, cd, thrust, air_density = aero_dict['max_q']['velocity'], aero_dict['max_q']['cd'], aero_dict['max_q']['thrust'], aero_dict['max_q']['air_density']
    dragNose = 0.5 * cd * air_density * S * velocity**2
    
    return

# Graph shear forces
def graphShear(shear_array, totalLength):
    dx = totalLength / len(shear_array)
    x = [i * dx for i in range(len(shear_array))]
    plt.plot(x, shear_array)
    plt.xlabel("Distance from aft (m)")
    plt.ylabel("Shear Force (N)")
    plt.title("Shear Forces")
    plt.show()

# Graph bending forces
def graphBending(bending_array, totalLength):
    dx = totalLength / len(bending_array)
    x = [i * dx for i in range(len(bending_array))]
    plt.plot(x, bending_array)
    plt.xlabel("Distance from aft (m)")
    plt.ylabel("Bending Forces (N * m)")
    plt.title("Bending Forces")
    plt.show()

# sfd_inputs = pd.ExcelFile('sfd_inputs.xlsx', engine='openpyxl')
# rocket_dict = getRocketSections(sfd_inputs)
# point_masses = getPointMasses(sfd_inputs)
# aero_dict = getAeroProperties(sfd_inputs)
# mass_model = getMassModel(rocket_dict, point_masses, 0.005)
# cg = getCG(rocket_dict, point_masses, getTotalMass(rocket_dict))
# length = getTotalLength(rocket_dict)

# print(rocket_dict) # TEST
# print(point_masses) # TEST
'''
