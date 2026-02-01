'''
Pathfiinder 3DOF SFD and RDOF

Author: Gary Huang
Updated: February 1, 2026

Purpose: To create a 3DOF simulation that calculates internal forces (Shear, Bending, Axial)
         along the rocket for recovery (I'm remaking rdof)

Methodology: Begin with 

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
    return np.lingalg.norm(vector)

def normalize(vector):
    mag = magnitude(vector)
    if mag == 0:
        return vector
    return vector / mag

xhat = np.array([1, 0]) # Unit vector pointing in the direction of the rocket nose
yhat = np.array([0, 1]) # Unit vector perpendicular to the rocket nose
# ------------------------------------------------------------------------------

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

rocket_dict_dry = vehicle.rocket_dict_dry
parachute_mass = vehicle.parachute_mass  # [kg]
recovery_bay_start = rocket_dict_dry["recovery_bay"]["bottom_distance_from_aft"]  # [m]
drag_coefficient = 2.2 # [-]
canopy_area = ((14 * c.FT2M) / 2)**2 * np.pi # [m^2]

dt = 0.1 # [s] Time step
time_before_deploy = 2 # [s] Time before deployment

# Time before deployment
for t in np.arange(0, time_before_deploy, dt):
    velocity = velocity + gravity * dt # [m / s] Update velocity due to gravity
    position = position + velocity * dt # [m] Update position
    # print(f"Time: {t:.2f} s, Position: {position}, Velocity: {velocity}")

# Deployment
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

def inflation_time(canopy_factor, diamater, velocity):
    '''
    canopy_factor: Canopy factor (2.5 for low porosity canopies, will be used)
    diameter: Diameter of parachute [m]
    velocity: Velocity of rocket at deployment [m / s]
    inflation_time: Time for parachute to fully inflate [s]
    '''
    inflation_time = (canopy_factor * diamater) / (velocity**0.85)
    return inflation_time

linear_density_array, length_along_rocket_linspace = mass_model(rocket_dict_dry, parachute_mass)
dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]  # [m]
total_mass = np.sum(linear_density_array * dx) # [kg]
cg = sfd.calcCG(linear_density_array, length_along_rocket_linspace) # [m] Center of gravity from bottom of rocket
terminal_velocity = sfd.calcTerminalVelocity(total_mass, magnitude(gravity), drag_coefficient, air_density, canopy_area) # [m / s]

