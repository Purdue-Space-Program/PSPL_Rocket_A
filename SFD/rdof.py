import math
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import vehicle_parameters as vehicle
import sfd as sfd
import loads as loads # This will cause plots of loads.py to show up when rdof.py is run
import parseWind as pw
from scipy.io import savemat
# ------------------------------------------------------------------------------

# Constants
LBF2N = 4.44822  # Pounds force to Newtons
FT2M = 0.3048  # Feet to Meters
N2LBF = 0.224809 # Newtons to Pounds force
M2FT = 3.28084 # Meters to Feet
# ------------------------------------------------------------------------------

# Input Parameters 
rocket_dict_dry = vehicle.rocket_dict_dry
parachute_mass = vehicle.parachute_mass  # [kg]
recovery_bay_start = rocket_dict_dry["recovery_bay"]["bottom_distance_from_aft"]  # [m]
max_q_velocity = vehicle.parameters.six_DoF_max_velocity  # [m / s]
AOA = sfd.calcAOA(loads.max_q_wind_gust, vehicle.parameters.six_DoF_max_velocity) # [radians] # NEED
wind_gust_speed = pw.percentile_75_wind_speed # [m/s]
horizontal_velocity = max_q_velocity * np.sin(AOA) # [m / s] # NEED
gravity = 9.81 # [m / s^2]
air_density = 1.81 # [kg / m^3] at apogee # NEED
drag_coefficent = 2.2 # []
canopy_area = (14 * FT2M / 2)**2 * np.pi # [m^2]
max_height = vehicle.parameters.six_DoF_estimated_apogee  # [m]
# ------------------------------------------------------------------------------

# Mass Model
def mass_model(rocket_dict, parachute_mass):
    '''
    rocket_dict: Dictionary of rocket components
    parachute_mass: Mass of parachute [kg]
    linear_density_array: Array of linear density along the rocket [kg / m]
    length_along_rocket_linspace: Array of length along rocket [m]
    '''
    num_points = 500
    length_along_rocket_linspace = np.linspace(rocket_dict["engine"]["bottom_distance_from_aft"], rocket_dict["nosecone"]["bottom_distance_from_aft"] + rocket_dict["nosecone"]["length"], num_points)
    dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]  # [m]

    linear_density_array = np.zeros(num_points)

    for component in rocket_dict:
        mass = rocket_dict[component]["mass"]
        if component == "recovery_bay":
            mass -= parachute_mass
        bottom_distance_from_aft = rocket_dict[component]["bottom_distance_from_aft"]
        length = rocket_dict[component]["length"]
        linear_density = mass / length
        for index, length_along_rocket in enumerate(length_along_rocket_linspace):
            above_component_bottom = length_along_rocket >= bottom_distance_from_aft
            below_component_top = length_along_rocket <= (bottom_distance_from_aft + length)
            
            if (above_component_bottom and below_component_top):
                linear_density_array[index] += linear_density
    recovery_bay_start = rocket_dict["recovery_bay"]["bottom_distance_from_aft"]
    index = int(recovery_bay_start / dx)
    linear_density_array = linear_density_array[:index]
    length_along_rocket_linspace = length_along_rocket_linspace[:index]

    return linear_density_array, length_along_rocket_linspace

linear_density_array, length_along_rocket_linspace = mass_model(rocket_dict_dry, parachute_mass)
total_mass = np.sum(linear_density_array * (length_along_rocket_linspace[1] - length_along_rocket_linspace[0])) # [kg]
# ------------------------------------------------------------------------------

# Functions for calculations
# Calculate inertia
def calcRotationalInertia(linear_density_array, length_along_rocket_linspace, cg):
    '''
    linear_density_array: Array of linear density across rocket [array]
    length_along_rocket_linspace: Numpy linspace for lengths along rocket [array]
    cg: Location of center of gravity [m]
    inertia: Rotational inertia around center of gravity [kg m^2]
    '''
    dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]
    mass_model = linear_density_array * dx
    inertia = 0
    for x in range(len(length_along_rocket_linspace)):
        inertia += mass_model[x] * (length_along_rocket_linspace[x] - cg)**2
    return inertia

# Calculate center of gravity (nosecone removed)
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
    terminal_velocity = math.sqrt((2 * mass * gravity) / (cd * rho * area))
    return terminal_velocity

# Calculate lateral acceleration
def calcLateralAcceleration(drag_force, total_mass):
    '''
    drag_force: Parachute drag force [N]
    total_mass: Total mass of rocket [kg]
    ay: Lateral acceleration [m / s^2]
    '''
    ay = (drag_force / total_mass)
    return ay

# Calculate angular acceleration
def calcAngularAcceleration(drag_force, recovery_bay_start, inertia, cg):
    '''
    drag_force: Parachute drag force [N]
    recovery_bay_start: Location of start of recovery bay [m]
    inertia: Rotational inertia around center of gravity [kg m^2]
    cg: Location of center of gravity of rocket [m]
    r: Angular acceleration [radians / s^2]
    '''
    r = ((-1) * drag_force * (abs(recovery_bay_start - cg))) / inertia
    return r

# Calculate shear forces
def calcShear(drag_force, recovery_bay_start, ay, linear_density_array, length_along_rocket_linspace, r, cg):
    '''
    drag_force: Parachute drag force [N]
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
    shear_array[int(recovery_bay_start / dx) - 1:] += drag_force

    return shear_array

# Calculate bending moments
def calcBending(shear_array, length_along_rocket_linspace):
    '''
    shear_array: Array of shear forces across rocket length [N]
    length_along_rocket_linspace: Array of rocket lengths [m]
    bending_array: Array of bending forces across rocket length [N m]
    '''
    dy = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]
    bending_array = np.cumsum(shear_array) * dy
    return bending_array

# Calculate axial forces
def calcAxial(drag_force, linear_density_array, length_along_rocket_linspace):
    '''
    drag_force: Parachute drag force [N]
    linear_density_array: Array of linear density across rocket length [kg / m]
    length_along_rocket_linspace: Array of rocket lengths [m]
    axial_array: Array of axial forces across rocket length [N]
    '''
    dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]
    ax = drag_force / total_mass
    mass_model = np.cumsum(linear_density_array * dx) # aft to nose
    axial_array = ax * mass_model

    return axial_array

# ------------------------------------------------------------------------------

# Calculations
terminal_velocity = calcTerminalVelocity(total_mass, gravity, drag_coefficent, air_density, canopy_area) # [m/s] solving for velocity setting weight and Drag equal
drag_force = calcDragForce(drag_coefficent, air_density, horizontal_velocity + wind_gust_speed, canopy_area) # [N]
descent_time = max_height/terminal_velocity # [s]

cg = calcCG(linear_density_array, length_along_rocket_linspace) # [m] Center of gravity
inertia = calcRotationalInertia(linear_density_array, length_along_rocket_linspace, cg) # [kg m^2] Rotational inertia
ay = calcLateralAcceleration(drag_force, total_mass) # [m / s^2] Lateral acceleration
r = calcAngularAcceleration(drag_force, recovery_bay_start, inertia, cg) # [radians / s^2] Angular acceleration
shear_array = np.array(calcShear(drag_force, recovery_bay_start, ay, linear_density_array, length_along_rocket_linspace, r, cg)) # [N] Shear force array
bending_array = np.array(calcBending(shear_array, length_along_rocket_linspace)) # [N m] Bending moment array
axial_array = np.array(calcAxial(drag_force, linear_density_array, length_along_rocket_linspace)) # [N] Axial forces array
# ------------------------------------------------------------------------------

# Converting to matlab file
matlab_dict = {"axial_array": axial_array, "shear_array": shear_array, "bending_array": bending_array, "length_along_rocket_linspace": length_along_rocket_linspace} # Dictionary to save as .mat file
savemat("rfd_outputs_revcovery.mat", matlab_dict) # Save as .mat file for MATLAB
# ------------------------------------------------------------------------------

# Print outputs
print("Outputs at recovery:")
print(f"Max axial force at recovery: {max(axial_array) * N2LBF:.2f} lbf")
print(f"Max shear force at recovery: {max(shear_array) * N2LBF:.2f} lbf")
print(f"Max bending moment at recovery: {max(bending_array) * N2LBF * M2FT:.2f} lbf-ft")
print("-----------------------------------")

print("Inputs:")
print(f"Angle of attack from max_q: {AOA * (180 / np.pi):.2f} degrees")
print(f"Drag force: {drag_force:.2f} N")
print(f"Wind gust at recovery: {wind_gust_speed:.2f} m/s")
print(f"Horizontal velocity at recovery: {horizontal_velocity:.2f} m/s")
print(f"Total velocity used for drag force at recovery: {(horizontal_velocity + wind_gust_speed):.2f} m/s")
print(f"Air density at apogee: {air_density:.2f} kg/m^3")
print(f"Canopy area: {canopy_area:.2f} m^2")
print(f"Drag coefficient at recovery: {drag_coefficent:.2f}")
print(f"Rocket mass at recovery: {total_mass:.2f} kg")
print("-----------------------------------")

print("Calculated values:")
print(f"Terminal velocity: {terminal_velocity:.2f} m/s")
print(f"Descent time: {descent_time:.2f} seconds")
print("-----------------------------------")
# ------------------------------------------------------------------------------

# Plotting
plt.figure()
plot_num = 1
for variable in ["shear_array", "bending_array", "axial_array"]:
    if variable == "shear_array":
        plot = shear_array * N2LBF
        ylabel = "Shear Force [lbf]"
        title = f"Shear Forces at Recovery"
    if variable == "bending_array":
        plot = bending_array * N2LBF * M2FT
        ylabel = "Bending Moment [lbf-ft]"
        title = f"Bending Moments at Recovery"
    if variable == "axial_array":
        plot = axial_array * N2LBF
        ylabel = "Axial Force [lbf]"
        title = f"Axial Forces at Recovery"
    plt.subplot(1,3, plot_num)
    plt.plot(length_along_rocket_linspace * M2FT, plot)
    plt.title(title)
    plt.xlabel("Length from aft [ft]")
    plt.ylabel(ylabel)
    plt.grid()
    plot_num += 1
plt.show()
# ------------------------------------------------------------------------------


'''
plt.plot(length_along_rocket_linspace, bending_array)
plt.title("Bending Moment vs Length Along Rocket")
plt.xlabel("Length Along Rocket [m]")
plt.ylabel("Bending Moment [Nm]")
plt.grid()
plt.show()
'''