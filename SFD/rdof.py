import os
os.chdir(os.path.dirname(__file__))

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import vehicle_parameters as vehicle
import constants as c

try:
    from SFD import parseWind
    from SFD import sfd
    from SFD import loads # This will cause plots of loads.py to show up when rdof.py is run
except ModuleNotFoundError:
    import parseWind
    import sfd
    import loads # This will cause plots of loads.py to show up when rdof.py is run
# ------------------------------------------------------------------------------

def CalculateCircleArea(diameter):
    radius = diameter/2
    area = np.pi * (radius**2)
    return(area)

# Input Parameters 
rocket_dict_dry = vehicle.rocket_dict_dry
parachute_mass = vehicle.parachute_mass  # [kg]
recovery_bay_start = rocket_dict_dry["recovery_bay"]["bottom_distance_from_aft"]  # [m]
max_height = vehicle.parameters.six_DoF_estimated_apogee  # [m]

# Scuffed velocity
max_q_velocity = vehicle.parameters.six_DoF_max_velocity # [m/s]
AOA_max_q = sfd.calcAOA(loads.max_q_wind_gust, vehicle.parameters.six_DoF_max_velocity) # [radians] # NEED
wind_gust_speed = parseWind.percentile_75_wind_gust_speed # [m/s]
horizontal_velocity = max_q_velocity * np.sin(AOA_max_q) # [m/s] # NEED

# Constants
gravity = 9.81 # [m / s^2]
air_density = 1.225 # [kg / m^3]
drag_coefficient = 2.2 # [-]

# Parachute parameters
# canopy_area = ((14 * c.FT2M) / 2)**2 * np.pi # [m^2]
canopy_outer_diameter = 29.56 * c.FT2M # [m] Diameter of parachute for area calculation, Rocketman High Performance CD 2.2 Parachute: https://www.the-rocketman.com/products/rocketman-high-performance-cd-2-2-parachutes?variant=42195940114526
canopy_inner_diameter = 14 * c.FT2M # [m] Diameter of parachute for area calculation, Rocketman High Performance CD 2.2 Parachute: https://www.the-rocketman.com/products/rocketman-high-performance-cd-2-2-parachutes?variant=42195940114526
canopy_area = CalculateCircleArea(canopy_outer_diameter) - CalculateCircleArea(canopy_inner_diameter) # [m^2] 168 inch (diameter?) parachute: https://www.the-rocketman.com/products/rocketman-high-performance-cd-2-2-parachutes?variant=42195940114526
canopy_factor = 2.5 # [-] Canopy factor for low porosity canopies, will be used for inflation time calculation
# ------------------------------------------------------------------------------

# Calculate AOA_recovery
orientation_start = 0 # [radians] Starting orientation at apogee measured from the horizontal
time_before_deployment = 0 # [s] Time from apogee to parachute deployment
vertical_velocity_deploy = (-1) * gravity * time_before_deployment # [m / s] Vertical velocity at parachute deployment
velocity_deploy = np.sqrt(vertical_velocity_deploy**2 + (horizontal_velocity + wind_gust_speed)**2) # [m / s] Total relative velocity at parachute deployment
angle_from_horizontal_parachute_deploy = np.arctan(vertical_velocity_deploy / (horizontal_velocity + wind_gust_speed)) # [radians] Angle of rocket vs parachute velocity vector at parachute deployment
AOA_recovery = orientation_start - angle_from_horizontal_parachute_deploy # [radians] Angle of attack at parachute deployment
orientation_start_recovery_list = [AOA_recovery] # [radians] List of angles of attack at recovery for plotting
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

# Calculate lateral acceleration
def calcParachuteAcceleration(drag_force, total_mass):
    '''
    drag_force: Parachute drag force [N]
    total_mass: Total mass of rocket [kg]
    acceleration: Parachute induced acceleration [m / s^2]
    '''
    acceleration = (drag_force / total_mass)
    return acceleration

# Calculate angular acceleration
def calcAngularAcceleration(drag_force, recovery_bay_start, inertia, cg, angle):
    '''
    drag_force: Parachute drag force [N]
    recovery_bay_start: Location of start of recovery bay [m]
    inertia: Rotational inertia around center of gravity [kg m^2]
    cg: Location of center of gravity of rocket [m]
    angle: Angle of attack at recovery [radians]
    angular_acceleration: Angular acceleration [radians / s^2]
    '''
    angular_acceleration = ((-1) * (drag_force * np.sin(angle)) * (abs(recovery_bay_start - cg))) / inertia
    return angular_acceleration

# Calculate shear forces
def calcShear(drag_force, recovery_bay_start, ay, linear_density_array, length_along_rocket_linspace, r, cg, angle):
    '''
    drag_force: Parachute drag force [N]
    ay: Lateral acceleration [m / s^2]
    linear_density_array: Array of linear density across rocket length [kg / m]
    length_along_rocket_linspace: Array of rocket lengths [m]
    r: Angular acceleration [radians / s^2]
    cg: Location of center of gravity of rocket [m]
    angle: Angle of attack at recovery [radians]
    shear_array: Array of shear forces across rocket length [N]
    '''
    dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]
    cg_rel_lengths = np.array(-1 * length_along_rocket_linspace + cg)
    mass_model = np.cumsum(linear_density_array * dx) # aft to nose
    cumulative_moment_about_cg = np.cumsum(linear_density_array * dx * cg_rel_lengths) # aft to nose
    shear_array = (-1) * ay * np.sin(angle) * mass_model - r * cumulative_moment_about_cg
    shear_array[int(recovery_bay_start / dx) - 1:] += (drag_force * np.sin(angle))

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
def calcAxial(drag_force, linear_density_array, length_along_rocket_linspace, angle):
    '''
    drag_force: Parachute drag force [N]
    linear_density_array: Array of linear density across rocket length [kg / m]
    length_along_rocket_linspace: Array of rocket lengths [m]
    angle: Angle of attack at recovery [radians]
    axial_array: Array of axial forces across rocket length [N]
    '''
    dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]
    ax = (drag_force * np.cos(angle)) / total_mass
    mass_model = np.cumsum(linear_density_array * dx) # aft to nose
    axial_array = ax * mass_model

    return axial_array

# ------------------------------------------------------------------------------

# Calculations
terminal_velocity = calcTerminalVelocity(total_mass, gravity, drag_coefficient, air_density, canopy_area) # [m/s] solving for velocity setting weight and Drag equal
drag_force = calcDragForce(drag_coefficient, air_density, velocity_deploy, canopy_area) # [N]
descent_time = max_height/terminal_velocity # [s]

inflation_time = sfd.inflation_time(canopy_factor, canopy_outer_diameter, velocity_deploy) # [s]
open_rate = canopy_area / inflation_time # [m^2 / s] Rate at which canopy area opens
dt = 0.1 # [s]


cg = calcCG(linear_density_array, length_along_rocket_linspace) # [m] Center of gravity
parachute_to_cg = recovery_bay_start  - cg # [m] Distance from parachute attachment point to CG
inertia = calcRotationalInertia(linear_density_array, length_along_rocket_linspace, cg) # [kg m^2] Rotational inertia
ay = calcParachuteAcceleration(drag_force, total_mass) # [m / s^2] Lateral acceleration

worst_shear_angle = np.pi / 2
worst_axial_angle = 0
r_worst_shear = calcAngularAcceleration(drag_force, recovery_bay_start, inertia, cg, worst_shear_angle) # [radians / s^2] Angular acceleration for worst shear angle
r_worst_axial = calcAngularAcceleration(drag_force, recovery_bay_start, inertia, cg, worst_axial_angle) # [radians / s^2] Angular acceleration for worst axial angle
worst_shear_array = np.array(calcShear(drag_force, recovery_bay_start, ay, linear_density_array, length_along_rocket_linspace, r_worst_shear, cg, worst_shear_angle)) # [N] Shear force array
worst_bending_array = np.array(calcBending(worst_shear_array, length_along_rocket_linspace)) # [N m] Bending moment array
worst_axial_array = np.array(calcAxial(drag_force, linear_density_array, length_along_rocket_linspace, worst_axial_angle)) # [N] Axial forces array
# ------------------------------------------------------------------------------

# Converting to matlab file
# matlab_dict = {"axial_array": axial_array, "shear_array": shear_array, "bending_array": bending_array, "length_along_rocket_linspace": length_along_rocket_linspace} # Dictionary to save as .mat file
# savemat("rfd_outputs_recovery.mat", matlab_dict) # Save as .mat file for MATLAB
# ------------------------------------------------------------------------------

if __name__ == "__main__":
    # Print outputs
    print("Outputs at recovery:")
    print(f"Worst axial force at recovery (0 degrees AOA): {max(worst_axial_array) * c.N2LBF:.2f} lbf")
    print(f"Worst shear force at recovery (90 degrees AOA): {max(worst_shear_array) * c.N2LBF:.2f} lbf")
    print(f"Worst bending moment at recovery (90 degrees AOA): {max(worst_bending_array) * c.N2LBF * c.M2FT:.2f} lbf-ft")
    print("-----------------------------------")

    print("Inputs:")
    print(f"Angle of attack from max_q: {AOA_max_q * (180 / np.pi):.2f} degrees")
    print(f"Orientation measured from horizontal at recovery: {orientation_start * (180 / np.pi):.2f} degrees")
    print(f"Angle of attack at recovery: {AOA_recovery * (180 / np.pi):.2f} degrees")
    print(f"Angle of parachute from horizontal: {angle_from_horizontal_parachute_deploy * (180 / np.pi):.2f} degrees")
    print(f"Drag force: {drag_force:.2f} N")
    print(f"Wind gust at recovery: {wind_gust_speed:.2f} m/s")
    print(f"Horizontal velocity at recovery: {horizontal_velocity:.2f} m/s")
    print(f"Vertical velocity at recovery: {vertical_velocity_deploy:.2f} m/s")
    print(f"Total velocity used for drag force at recovery: {velocity_deploy:.2f} m/s")
    print(f"Air density at apogee: {air_density:.2f} kg/m^3")
    print(f"Canopy area: {canopy_area:.2f} m^2")
    print(f"Drag coefficient at recovery: {drag_coefficient:.2f}")
    print(f"Rocket mass at recovery: {total_mass:.2f} kg")
    print("-----------------------------------")

    print("Calculated values:")
    print(f"Terminal velocity: {terminal_velocity:.2f} m/s")
    print(f"Descent time: {descent_time:.2f} seconds")
    print("-----------------------------------")
# ------------------------------------------------------------------------------

# Plotting and converting to matlab file
matlab_dict = {}
plt.figure()
plot_num = 1

limit_load_shear_array = None
limit_load_bending_array = None
limit_load_axial_array = None

for orientation in orientation_start_recovery_list:
    angle = orientation - angle_from_horizontal_parachute_deploy
    for t in np.arange(0, inflation_time, dt):
        area = open_rate * t
        if area > canopy_area:
            area = canopy_area
        drag_force = calcDragForce(drag_coefficient, air_density, velocity_deploy, area) # [N]
        lateral_acceleration = calcParachuteAcceleration(drag_force, total_mass) # [m / s^2]
        # Orientation update at the end
        # Velocity update at the end



    r_i = calcAngularAcceleration(drag_force, recovery_bay_start, inertia, cg, angle) # [radians / s^2] Angular acceleration if rocket at angle at recovery
    shear_array = np.array(calcShear(drag_force, recovery_bay_start, ay, linear_density_array, length_along_rocket_linspace, r_i, cg, angle)) # [N] Shear force array if rocket at angle at recovery
    bending_array = np.array(calcBending(shear_array, length_along_rocket_linspace)) # [N m] Bending moment array if rocket at angle at recovery
    axial_array = np.array(calcAxial(drag_force, linear_density_array, length_along_rocket_linspace, angle)) # [N] Axial forces array if rocket at angle at recovery
    
    if (limit_load_shear_array is None) or (np.mean(shear_array) > np.mean(limit_load_shear_array)):
        limit_load_shear_array = shear_array
    
    if (limit_load_bending_array is None) or (np.mean(bending_array) > np.mean(limit_load_bending_array)):
        limit_load_bending_array = bending_array
    
    if (limit_load_axial_array is None) or (np.mean(axial_array) > np.mean(limit_load_axial_array)):
        limit_load_axial_array = axial_array
    # matlab_dict[f"AOA_recovery_deg"] = angle * (180 / np.pi)
    
    plot_def = [
        (shear_array, c.N2LBF, "Shear Force [lbf]", "Shear Forces"),
        (bending_array, c.N2LBF * c.M2FT, "Bending Moment [lbf-ft]", "Bending Moments"),
        (axial_array, c.N2LBF, "Axial Force [lbf]", "Axial Forces")
    ]

    for data, scale, ylabel, base_title in plot_def:
        plot = data * scale
        title = f"{base_title} at {angle * (180 / np.pi):.0f} deg at recovery"
        plt.subplot(len(orientation_start_recovery_list), 3, plot_num)
        plt.plot(length_along_rocket_linspace * c.M2FT, plot)
        plt.title(title)
        plt.xlabel("Length from aft [ft]")
        plt.ylabel(ylabel)
        
        plt.ylim(top = 1.2 * max(max(plot), 100), bottom = min(-abs(min(plot)) * 1.2, -100)) # to show low values when values are like 10^-13 
        
        plt.grid()
        plot_num += 1


matlab_dict[f"shear_array_recovery"] = limit_load_shear_array
matlab_dict[f"bending_array_recovery"] = limit_load_bending_array
matlab_dict[f"axial_array_recovery"] = limit_load_axial_array
matlab_dict["length_along_rocket_linspace"] = length_along_rocket_linspace

savemat("rfd_outputs_recovery.mat", matlab_dict) # Save as .mat file for MATLAB
# ------------------------------------------------------------------------------

if len(orientation_start_recovery_list) > 1:
    plt.tight_layout()
# ------------------------------------------------------------------------------

if __name__ == "__main__":
    plt.show()
print(drag_force)