'''
Pathfiinder 3DOF SFD and RDOF

Author: Gary Huang
Updated: February 1, 2026

Purpose: To create a 3DOF simulation that calculates internal forces (Shear, Bending, Axial)
         along the rocket for recovery (I'm remaking rdof)

Methodology: Get attitude of rocket at each time step for AOA, calculate drag force, calculate net force, update velocity and position, 
calculate torque, update orientation, repeat until landing. Then use the time history of forces and orientations to calculate internal forces along the rocket.

Inputs: 
    Initial position, velocity, thrust force, initial linear density array, initial length along rocket linspace,
    air density, drag coefficient, cross sectional area, wind gust speed
Outputs:
    Shear, Bending, Axial force arrays saved as .mat files for MATLAB
'''

import numpy as np

import matplotlib.pyplot as plt
import sfd as sfd
import constants as c
from scipy.io import savemat
from scipy.integrate import solve_ivp
import parseWind as pw

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import vehicle_parameters as vehicle
# ------------------------------------------------------------------------------

# Vector operations
def magnitude(vector):
    return np.linalg.norm(vector)

def normalize(vector):
    mag = magnitude(vector)
    if mag == 0:
        return vector
    return vector / mag

xhat = np.array([1, 0]) # Unit vector pointing in the direction of the rocket nose
yhat = np.array([0, 1]) # Unit vector perpendicular to the rocket nose
# ------------------------------------------------------------------------------

# Simulation Storage
time_array = []
position_array = []
velocity_array = []
acceleration_array = []
orientation_array = []
AOA_array = []

torque_array = []
angular_acceleration_array = []
angular_velocity_array = []

F_drag_array = [] # [N] Array to store drag forces for each time step
area_array = [] # [m^2] Array to store canopy area for each time step
lateral_acceleration_array = [] # [m / s^2] Array to store lateral acceleration for each time step

shear_force_array = [] # [N] Array to store shear forces along the rocket for each time step
bending_moment_array = [] # [N * m] Array to store bending moments along the rocket for each time step
axial_force_array = [] # [N] Array to store axial forces along the rocket for each time step

# Initial conditions
gravity = 9.81 * -yhat # [m / s^2]
air_density = 1.225 # [kg / m^3]
AOA = 0 # [radians]

# VERY SCUFFED VELOCITY 
max_q_velocity = vehicle.parameters.six_DoF_max_velocity # [m / s] Velocity at max q
wind_gust = pw.percentile_75_wind_gust_speed # [m / s] Assume constant wind gust
max_q_AOA = sfd.calcAOA(wind_gust, max_q_velocity) # [radians] Angle of attack at max q
velocity = max_q_velocity * np.sin(max_q_AOA) * xhat # [m / s] Adjusted velocity for AOA
max_height = vehicle.parameters.six_DoF_estimated_apogee  # [m]
position = np.array([0, max_height])  # [m] Initial position at apogee

orientation = [1, 0] # Initial orientation of rocket in polar coordinates [magnitude, angle]

rocket_dict_dry = vehicle.rocket_dict_dry # Dictionary of rocket components when dry
parachute_mass = vehicle.parachute_mass  # [kg]
recovery_bay_start = rocket_dict_dry["recovery_bay"]["bottom_distance_from_aft"]  # [m]
drag_coefficient = 2.2 # [-]
canopy_area = ((14 * c.FT2M) / 2)**2 * np.pi # [m^2]

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
            below_component_top = length_along_rocket < (bottom_distance_from_aft + length)
            
            if (above_component_bottom and below_component_top):
                linear_density_array[index] += linear_density
    recovery_bay_start = rocket_dict["recovery_bay"]["bottom_distance_from_aft"]
    index = int(recovery_bay_start / dx)
    linear_density_array = linear_density_array[:index]
    length_along_rocket_linspace = length_along_rocket_linspace[:index]

    return linear_density_array, length_along_rocket_linspace

linear_density_array, length_along_rocket_linspace = mass_model(rocket_dict_dry, parachute_mass)
dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]  # [m]
total_mass = np.sum(linear_density_array * dx) # [kg]
cg = sfd.calcCG(linear_density_array, length_along_rocket_linspace) # [m] Center of gravity from bottom of rocket
parachute_to_cg = recovery_bay_start  - cg # [m] Distance from parachute attachment point to CG
terminal_velocity = sfd.calcTerminalVelocity(total_mass, magnitude(gravity), drag_coefficient, air_density, canopy_area) # [m / s]
inertia = sfd.calcRotationalInertia(linear_density_array, length_along_rocket_linspace, cg) # [kg * m^2] Rotational inertia about CG
angular_velocity = 0 # [rad / s] Initial angular velocity
F_g = total_mass * gravity # [N] Gravitational force

inflation_duration = sfd.inflation_time(2.5, np.sqrt(canopy_area / np.pi) * 2, magnitude(velocity)) # [s]
open_rate = canopy_area / inflation_duration # [m^2 / s] Rate at which canopy area opens

dt = 0.1 # [s] Time step
time_before_deploy = 2 # [s] Time before deployment

# Time before deployment
for t in np.arange(0, time_before_deploy, dt):
    position_array.append(position)
    velocity_array.append(velocity)    
    velocity = velocity + gravity * dt # [m / s] Update velocity due to gravity
    position = position + velocity * dt # [m] Update position
    time_array.append(t)
    acceleration_array.append(gravity)
    orientation_array.append(orientation)
    AOA_array.append(AOA)

    torque_array.append(0) # No torque before deployment
    angular_acceleration_array.append(0) # No angular acceleration before deployment
    angular_velocity_array.append(0) # No angular velocity before deployment
    
    F_drag_array.append(np.array([0, 0])) # No drag force before deployment
    area_array.append(0) # No canopy area before deployment
    lateral_acceleration_array.append(0) # No lateral acceleration before deployment
    shear_force_array.append(np.zeros_like(linear_density_array)) # No shear force before deployment
    bending_moment_array.append(np.zeros_like(linear_density_array)) # No bending moment before deployment\
    axial_force_array.append(np.zeros_like(linear_density_array)) # No axial force before deployment
    
    # print(f"Time: {t:.2f} s, Position: {position}, Velocity: {velocity}")

# Deploying to total inflation
for t in np.arange(0, inflation_duration, dt):
    # Angle of attack update
    orientation_cartesian = orientation[0] * np.array([np.cos(orientation[1]), np.sin(orientation[1])]) # Convert orientation to cartesian
    AOA = np.arccos(np.dot(normalize(orientation_cartesian), normalize(velocity))) # [radians] Angle of attack between velocity vector and rocket orientation

    # Canopy area update
    area = open_rate * t # [m^2] Canopy area at time t

    # Force update
    F_drag = sfd.calcDragForce(drag_coefficient, air_density, magnitude(velocity), area) * (-1) * normalize(velocity) # [N] Drag force
    F_net = F_drag + F_g # [N] Net force

    # Acceleration Velocity Position update
    acceleration = F_net / total_mass # [m / s^2] Acceleration
    velocity = velocity + acceleration * dt # [m / s] Update velocity
    position = position + velocity * dt # [m] Update position

    # Torque update
    torque = np.cross(parachute_to_cg * orientation_cartesian, F_drag) # [N * m] Torque about CG
    angular_acceleration = magnitude(torque) / inertia # [rad / s^2] Angular acceleration
    angular_velocity = angular_acceleration * dt # [rad / s] Angular velocity
    orientation = [orientation[0], orientation[1] + angular_velocity * dt] # Update orientation angle

    # Calculate internal forces
    # Shear forces
    lateral_acceleration = np.abs(magnitude(F_drag) * np.sin(AOA) / total_mass) # [m / s^2] Lateral acceleration
    lateral_acceleration_array.append(lateral_acceleration) # Store lateral acceleration for this time step
    F_drag_perpendicular = magnitude(F_drag) * np.sin(AOA) # [N] Component of drag force perpendicular to rocket orientation
    
    shear_force = sfd.calcShearRecovery(F_drag_perpendicular, recovery_bay_start, lateral_acceleration, linear_density_array, length_along_rocket_linspace, angular_acceleration, cg)

    # Bending moments
    bending_moment = sfd.calcBending(shear_force, length_along_rocket_linspace)

    # Axial forces
    F_drag_parallel = magnitude(F_drag) * np.cos(AOA) # [N] Component of drag force parallel to rocket orientation
    axial_force = sfd.calcAxialRecovery(F_drag_parallel, linear_density_array, length_along_rocket_linspace, total_mass)

    '''
    F_drag = sfd.calcDragForce(drag_coefficient, air_density, magnitude(velocity), canopy_area) * (-1) * normalize(velocity) # [N] Drag force
    F_net = F_drag + F_g # [N] Net force
    
    acceleration = F_net / total_mass # [m / s^2] Acceleration
    velocity = velocity + acceleration * dt # [m / s] Update velocity
    position = position + velocity * dt # [m] Update position
    
    orientation_cartesian = orientation[0] * np.array([np.cos(orientation[1]), np.sin(orientation[1])]) # Convert orientation to cartesian
    
    

    torque = np.cross(parachute_to_cg * orientation_cartesian, magnitude(F_drag) * normalize(velocity)) # [N * m] Torque about CG
    angular_acceleration = magnitude(torque) / inertia # [rad / s^2] Angular acceleration
    angular_velocity = angular_acceleration * dt # [rad / s] Angular velocity
    AOA = AOA + angular_velocity * dt # [radians] Update AOA
    
    orientation = [orientation[0], orientation[1] + angular_velocity * dt] # Update orientation
    '''
    time_array.append(time_array[-1] + dt)
    position_array.append(position)
    velocity_array.append(velocity)
    acceleration_array.append(acceleration)
    orientation_array.append(orientation)
    AOA_array.append(AOA)

    torque_array.append(torque)
    angular_acceleration_array.append(angular_acceleration)
    angular_velocity_array.append(angular_velocity)

    area_array.append(area)
    F_drag_array.append(F_drag) # Store drag force for this time step
    shear_force_array.append(shear_force) # Store shear force array for this time step
    bending_moment_array.append(bending_moment) # Store bending moment array for this time step
    axial_force_array.append(axial_force) # Store axial force array for this time step


max_shear_forces = []
max_bending_moments = []
max_axial_forces = []
for i in range(len(time_array)):
    max_shear_forces.append(max(np.abs(shear_force_array[i])))
    max_bending_moments.append(max(np.abs(bending_moment_array[i])))
    max_axial_forces.append(max(np.abs(axial_force_array[i])))
worst_time_shear = max_shear_forces.index(max(max_shear_forces))
worst_time_bending = max_bending_moments.index(max(max_bending_moments))
worst_time_axial = max_axial_forces.index(max(max_axial_forces))





table = ""
for i in range(len(time_array)):
    table += f"Time: {time_array[i]:.2f} s, Position: {position_array[i]}, Velocity: {velocity_array[i]}, Acceleration: {acceleration_array[i]}, Orientation: {orientation_array[i]}, AOA: {AOA_array[i]}, Drag Force: {F_drag_array[i] if i < len(F_drag_array) else 0}, Canopy Area: {area_array[i] if i < len(area_array) else 0}\n"
print(table)

#print("Shear Force Array:")
#for i, sf in enumerate(shear_force_array):
#    print(f"Time: {time_array[i]:.2f} s, Shear Force: {sf}")
#print("Bending Moment Array:")
#for i, bm in enumerate(bending_moment_array):
#    print(f"Time: {time_array[i]:.2f} s, Bending Moment: {bm}")
#print("Axial Force Array:")
#for i, af in enumerate(axial_force_array):
#    print(f"Time: {time_array[i]:.2f} s, Axial Force: {af}")
'''
print("Time")
print((time_array))
print("Position")
print((position_array))
print("Velocity")
print((velocity_array))
print("Acceleration")
print((acceleration_array))
print("Orientation")
print((orientation_array))
print("AOA")
print((AOA_array))
print("Drag Force")
print((F_drag_array))'''
#print(max_shear_forces)
#print(len(max_bending_moments))
#print(len(max_axial_forces))
#print(f"Worst shear index: {worst_time_shear}, Worst bending index: {worst_time_bending}, Worst axial index: {worst_time_axial}")

# Plotting
worst_shear_array = shear_force_array[worst_time_shear]
worst_bending_array = bending_moment_array[worst_time_bending]
worst_axial_array = axial_force_array[worst_time_axial]
worst_cases = [worst_shear_array, worst_bending_array, worst_axial_array]

plot_on = False
if plot_on == True:
    plt.figure()
    plot_num = 1
    for variable in worst_cases:
        if variable is worst_shear_array:
            plot = variable * c.N2LBF
            ylabel = "Shear Force [lbf]"
            title = f"Worst Shear Forces at {time_array[worst_time_shear]:.2f} s"
        if variable is worst_bending_array:
            plot = variable * c.N2LBF * c.M2FT
            ylabel = "Bending Moment [lbf-ft]"
            title = f"Worst Bending Moments at {time_array[worst_time_bending]:.2f} s"
        if variable is worst_axial_array:
            plot = variable * c.N2LBF
            ylabel = "Axial Force [lbf]"
            title = f"Worst Axial Forces at {time_array[worst_time_axial]:.2f} s"
        plt.subplot(1,3, plot_num)
        plt.plot(length_along_rocket_linspace * c.M2FT, plot)
        plt.title(title)
        plt.xlabel("Length from aft [ft]")
        plt.ylabel(ylabel)
        plt.grid()
        plot_num += 1
    if __name__ == "__main__":
        plt.show()

print(magnitude(F_drag_array[22]) * np.sin(AOA_array[22]))