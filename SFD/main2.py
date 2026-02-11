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

# Vector operations
def magnitude(vector):
    sum = 0
    for component in vector:
        sum += component**2
    sum = np.sqrt(sum)
    return sum

def normalize(vector):
    mag = magnitude(vector)
    if mag == 0:
        return vector
    return vector / mag

rocket_dict_dry = vehicle.rocket_dict_dry
parachute_mass = vehicle.parachute_mass  # [kg]
recovery_bay_start = rocket_dict_dry["recovery_bay"]["bottom_distance_from_aft"]  # [m]
max_height = vehicle.parameters.six_DoF_estimated_apogee  # [m]

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
dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]  # [m]
total_mass = np.sum(linear_density_array * (length_along_rocket_linspace[1] - length_along_rocket_linspace[0])) # [kg]
cg = sfd.calcCG(linear_density_array, length_along_rocket_linspace) # [m] Center of gravity
inertia = sfd.calcRotationalInertia(linear_density_array, length_along_rocket_linspace, cg) # [kg*m^2] Rotational inertia
# ------------------------------------------------------------------------------

# Import values
apogee = vehicle.parameters.six_DoF_estimated_apogee
max_q_velocity = vehicle.parameters.six_DoF_max_velocity
wind_gust = pw.percentile_75_wind_gust_speed
max_AOA = np.arctan(wind_gust / max_q_velocity)  # Maximum angle of attack based on wind gust and max velocity

# Canopy parameters
canopy_outer_diameter = 14 * c.FT2M  # Convert canopy diameter from feet to meters
canopy_inner_diameter = (29.56 / 12) * c.FT2M  # Convert canopy diameter from inches to meters
canopy_area = np.pi * (canopy_outer_diameter / 2)**2 - np.pi * (canopy_inner_diameter / 2)**2  # Area of the canopy [m^2]
canopy_factor = 2.5 # [-] Canopy factor for low porosity canopies, will be used for inflation time calculation

air_density = 1.225  # Air density at sea level [kg/m^3]
drag_coefficient = 2.2  # Drag coefficient for a parachute, can vary
gravity = np.array([0.0, -9.81, 0.0])  # Gravity vector [m/s^2]

# Initial conditions
initial_position = np.array([0.0, vehicle.parameters.six_DoF_estimated_apogee, 0.0])  # Starting at the estimated apogee [m]
initial_velocity = np.array([max_q_velocity * np.sin(max_AOA), 0.0, 0.0])  # Initial velocity vector based on max AOA and max velocity [m/s]
initial_acceleration = gravity  # Initial acceleration (gravity) [m/s^2]
initial_angular_acceleration = 0.0  # Initial angular acceleration [rad/s^2]
initial_angular_velocity = 0.0  # Initial angular velocity [rad/s]
initial_orientation = 0.0 # Initial orientation [radians]

print(f"Initial Position: {initial_position}")
print(f"Initial Velocity: {initial_velocity}")
print(f"Initial Acceleration: {initial_acceleration}")
print(f"Initial Angular Acceleration: {initial_angular_acceleration} rad/s^2")
print(f"Initial Angular Velocity: {initial_angular_velocity} rad/s")
print(f"Initial Orientation: {initial_orientation} radians")
print(f"Total mass: {total_mass:.2f} kg")
print("--------------------------------------------------\n")

# Set up variables
position = initial_position
velocity = initial_velocity
acceleration = initial_acceleration
angular_acceleration = initial_angular_acceleration
angular_velocity = initial_angular_velocity
orientation = initial_orientation
time_dict = {}

# Time stepping before parachute deployment
time_before_deploy = 0
time = 0
dt = 0.1
time_steps_1 = int(time_before_deploy / dt)

for t in range(time_steps_1):
    print(f"Time: {time:.2f} s\n Position: {position}\n Velocity: {velocity}\n Acceleration: {acceleration}\n Orientation: {orientation:.2f} radians\n-----------------------------\n")
    # Update position and velocity using kinematic equations
    position += velocity * dt + 0.5 * acceleration * dt**2
    velocity += acceleration * dt
    time += dt

inflation_time = canopy_factor * canopy_outer_diameter / magnitude(velocity)**(0.85)
open_rate = canopy_area / inflation_time
area = 0

print(f"Inflation Time: {inflation_time:.4f} s")
print(f"Open Rate: {open_rate:.2f} m^2/s")
print(f"Area: {area:.2f} m^2")

time_steps_2 = int(inflation_time / dt)
time_array = []
for t in range(time_steps_2):
    array_i = [] # Time, velocity, shear array, axial array
    array_i.append(time)
    array_i.append(velocity)
    print(f"Time: {(time):.2f} s\n Position: {position}\n Velocity: {velocity}\n Acceleration: {acceleration}\n Orientation: {orientation:.2f} radians\n Area: {area:.2f} m^2\n")
    print(f"Angular_acceleration: {angular_acceleration:.4f} rad/s^2\n Angular_velocity: {angular_velocity:.4f} rad/s\n")
    
    F_drag = 0.5 * air_density * magnitude(velocity)**2 * drag_coefficient * area * normalize(velocity) * (-1)  # Drag force vector [N]
    
    # Update torque and angular acceleration based on drag force
    orientation_cartesian = np.array([np.cos(orientation), np.sin(orientation), 0.0])  # Orientation vector in Cartesian coordinates
    orientation_perpendicular = np.array([np.cos(orientation + np.pi/2), np.sin(orientation + np.pi/2), 0.0])  # Perpendicular vector to orientation

    print(f"Orientation Cartesian: {orientation_cartesian}\nOrientation Perpendicular: {orientation_perpendicular}\n")

    F_drag_perpendicular = np.dot(F_drag, orientation_perpendicular) * orientation_perpendicular  # Component of drag force perpendicular to orientation
    F_drag_parallel = np.dot(F_drag, orientation_cartesian) * orientation_cartesian  # Component of drag force parallel to orientation
    print(f"F_drag_perpendicular: {F_drag_perpendicular}\nF_drag_parallel: {F_drag_parallel}\n")
    
    lateral_acceleration = F_drag_perpendicular / total_mass  # Lateral acceleration due to perpendicular drag force [m/s^2]
    parallel_acceleration = F_drag_parallel / total_mass  # Parallel acceleration due to parallel drag force [m/s^2]
    print(f"Lateral Acceleration: {lateral_acceleration}\nParallel Acceleration: {parallel_acceleration}\n")
    
    # Calculate internal forces
    cg_relative_lengths = np.array(-1 * length_along_rocket_linspace + cg)  # Length from each segment to center of gravity [m]
    relative_lengths_index = 0
    moments_array = []
    mass_array = linear_density_array * (length_along_rocket_linspace[1] - length_along_rocket_linspace[0])  # Mass of each segment [kg]
    cumsum_mass_array = []
    slice_sum = 0
    for slice in mass_array:
        slice_sum += slice
        moment = slice_sum * cg_relative_lengths[relative_lengths_index]  # Moment contribution from each segment [N*m]
        moments_array.append(moment)
        cumsum_mass_array.append(slice_sum) # Validated against np.cumsum(mass_array)
        relative_lengths_index += 1
    print(f"Slice sum: {slice_sum}\n")
    #print(f"Moments Array: {moments_array}\n")
    #print(f"Validate moments array: {np.cumsum(linear_density_array * (length_along_rocket_linspace[1] - length_along_rocket_linspace[0]) * cg_relative_lengths)}\n")
    moments_array = np.array(moments_array)
    cumsum_mass_array = np.array(cumsum_mass_array)
    shear_array_i = magnitude(lateral_acceleration) * cumsum_mass_array + angular_acceleration * moments_array  # Shear force at each segment [N]
    shear_array_i[int(recovery_bay_start / dx) - 1:] += magnitude(F_drag_perpendicular)
    array_i.append(shear_array_i)
    
    axial_array_i = magnitude(parallel_acceleration) * cumsum_mass_array  # Axial force at each segment [N]
    array_i.append(axial_array_i)
    time_array.append(array_i)

    # Update states
    torque = magnitude(F_drag_perpendicular) * (recovery_bay_start - cg)  # Torque due to drag force [N*m]
    # Update angular acceleration based on torque
    angular_acceleration = torque / inertia  # [rad/s^2]
    # Update angular velocity and orientation
    angular_velocity += angular_acceleration * dt
    orientation += angular_velocity * dt

    # Update acceleration based on drag force
    acceleration = gravity + F_drag / total_mass

    # Update position and velocity using kinematic equations
    position += velocity * dt + 0.5 * acceleration * dt**2
    velocity += acceleration * dt
    time += dt

    print(f"Drag Force: {F_drag} N\n")
    print(f"Torque: {torque} N*m\n-----------------------------\n")
    area += open_rate * dt
    if area > canopy_area:
        area = canopy_area


plt.figure()
plt.plot(length_along_rocket_linspace, time_array[3][3], label="Axial Force")
plt.grid()
plt.show()

# Max axial force at t = 0.3 time_array[3][3]