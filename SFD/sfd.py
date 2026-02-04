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

# Make linear_dnsity_array and length_along_rocket_linspace inputs
def mass_model(rocket_dict):
    '''
    rocket_dict: Dictionary of rocket components
    linear_density_array: Array of linear density along the rocket [kg / m]
    length_along_rocket_linspace: Array of length along rocket [m]
    '''

    num_points = 500
    length_along_rocket_linspace = np.linspace(rocket_dict["engine"]["bottom_distance_from_aft"], rocket_dict["nosecone"]["bottom_distance_from_aft"] + rocket_dict["nosecone"]["length"], num_points)
    linear_density_array = np.zeros(num_points)
    for component in rocket_dict:
        mass = rocket_dict[component]["mass"]
        bottom_distance_from_aft = rocket_dict[component]["bottom_distance_from_aft"]
        length = rocket_dict[component]["length"]
        linear_density = mass / length
        for index, length_along_rocket in enumerate(length_along_rocket_linspace):
            above_component_bottom = length_along_rocket >= bottom_distance_from_aft
            below_component_top = length_along_rocket <= (bottom_distance_from_aft + length)
            
            if (above_component_bottom and below_component_top):
                linear_density_array[index] += linear_density
    return linear_density_array, length_along_rocket_linspace

# Calculate dynamic pressure
def calcQ(rho, velocity):
    '''
    rho: air density [kg / m^3]
    velocity: air velocity = rocket velocity [m/s]
    Q: Dynamic pressure [Pa]
    '''
    return rho * velocity**2 / 2

# Calculate angle of attack
def calcAOA(wind_gust, velocity):
    '''
    wind_gust: wind gust velocity [m/s]
    velocity: rocket velocity [m/s]
    AOA: Angle of attack [radians]
    '''
    # print(f"AOA: {degrees(atan(wind_gust / velocity))} degrees")
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
    # print(f"Fin mid-chord: {mid_chord}") # TEST
    Kfb = 1 + (diameter / 2) / (fin_height + diameter / 2) # Coefficient, Cambridge equation 32
    # print(f"Fin Kfb: {Kfb}") # TEST
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
def calcNoseSD(cg, noseCP, diameter):
    '''
    cg: Location of center of gravity of rocket [m]
    noseCP: Location of center of pressure of the nose [m]
    diameter: Rocket diameter [m]
    noseSD: Nose stability derivative [unitless]
    '''
    noseSD = abs(cg - noseCP) / diameter 
    return noseSD

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
    #mass_flow_rate = vehicle.parameters.mass_flow_rate
    #OF_ratio = vehicle.parameters.OF_ratio
    ox_flow_rate = vehicle.parameters.oxidizer_mass_flow_rate #mass_flow_rate * OF_ratio / (1 + OF_ratio)
    fuel_flow_rate = vehicle.parameters.fuel_mass_flow_rate #mass_flow_rate * 1 / (1 + OF_ratio)
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
            
        # print(f"{float(t)}, {float(vehicle_length - cg[-1])}")
    
    cg_max_q = cg[-1]
    x = int(0.4 / dt) # Hardcoded off the rail time at 0.4s
    cg_off_the_rail = cg[x]
    return cg, cg_max_q, cg_off_the_rail

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
    noseconeToFin: Length from nosecone to find [m]
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
    noseCP = 2/6 * nosecone_length # [m] Nosecone CP from nose tip
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
    masses = linear_density_array * dx # aft to nose
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

# Calculate drag force
def calcDragForce(cd, rho, velocity, area):
    '''
    cd: Drag coefficient
    rho: Air density [kg / m^3]
    velocity: Velocity [m / s]
    area: Reference area [m^2]
    drag_force: Parachute drag force [N]
    '''
    drag_force = 0.5 * cd * rho * (velocity ** 2) * area

    # print("*********************")
    # print(f"cd: {cd:.2f}")
    # print(f"rho: {rho:.2f}")
    # print(f"velocity: {velocity:.2f}")
    # print(f"area: {area:.2f}")

    return drag_force

# Calculate terminal velocity
def calcTerminalVelocity(mass, gravity, cd, rho, area):
    '''
    mass: Mass of rocket [kg]
    gravity: Gravitational acceleration [m / s^2]
    cd: Drag coefficient
    rho: Air density [kg / m^3]
    area: Reference area [m^2]
    terminal_velocity: Terminal velocity [m / s]
    '''
    terminal_velocity = np.sqrt((2 * mass * gravity) / (cd * rho * area))
    return terminal_velocity
