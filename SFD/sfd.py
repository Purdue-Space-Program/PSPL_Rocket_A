import pandas as pd
from math import ceil, sqrt, atan, pi, degrees
import numpy as np
import matplotlib.pyplot as plt


LB2KG = 0.453592
FT2M = 0.3048
IN2M = 0.0254
LBF2N = 4.44822
MPH2MPS = 0.44704

def codeFriendlyName(str, delimiter='_'):
    return delimiter.join(str.split(' (')[0].split(' ')).lower() # Delimiter is character between words ex. "_"

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
    return rho * velocity**2 / 2

# Calculate cross-sectional area
def getArea(diameter):
    return pi * (diameter / 2)**2

# Calculate lift force
def getLiftForce(Q, S, AOA, SD):
    return Q * S * AOA * SD # Aspire page 14

# Calculate rotational inertia
def getRotationalInertia(mass_model, cg, length):
    lengths = np.linspace(0, length, len(mass_model))
    rotationalInertia = np.sum(mass_model * (lengths - cg)**2)
    return rotationalInertia

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
    return sum / inertia

# Calculate shear force
def getShearForce(mass_model, totalLength, lift_dict, cp_dict, ay, r, cg):
    dy = totalLength / len(mass_model)
    lengths = np.linspace(0, totalLength, len(mass_model))
    # print(lengths) # TEST
    
    # Calculate shear force at each point on rocket Aspire page 21
    shear_array = (-1) * (ay * np.cumsum(mass_model) + r * np.cumsum(mass_model * (cg - lengths)))
    shear_array[int(cp_dict['nose'] / dy):] += lift_dict['nose']
    shear_array[int(cp_dict['fin'] / dy):] += lift_dict['fin']
    shear_array[int(cp_dict['boattail'] / dy):] -= lift_dict['boattail']
    
    # print(shear_array) # TEST
    return shear_array # [array]

# Calculate bending force
def getBendingForce(shear_array, totalLength):
    dy = totalLength / len(shear_array)
    bending_array = np.cumsum(shear_array) * dy # Bending graph is integral of shear graph
    return bending_array # [array]

def getAxialForces(aero_dict, mass_model):
    return

# Graph shear forces
def graphShear(shear_array, totalLength):
    dx = totalLength / len(shear_array)
    x = [i * dx for i in range(len(shear_array))]
    plt.plot(x, shear_array)
    plt.title("Shear Forces")
    plt.show()

# Graph bending forces
def graphBending(bending_array, totalLength):
    dx = totalLength / len(bending_array)
    x = [i * dx for i in range(len(bending_array))]
    plt.plot(x, bending_array)
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

