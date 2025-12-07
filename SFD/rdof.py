import math
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import vehicle_parameters as vehicle
import sfd as sfd
import loads as loads

linear_density_array, length_along_rocket_linspace = sfd.mass_model(vehicle.rocket_dict_recovery)

total_mass = np.sum(linear_density_array * (length_along_rocket_linspace[1] - length_along_rocket_linspace[0])) # [kg]
gravity = 9.81
air_density = 1.81
drag_coefficent = 2.2
canopy_area = 24.6677824063
max_height = vehicle.parameters.estimated_apogee  # meters

''' Sphereacutes specs '''
'''
rocket_mass = 69
gravity = 9.81
air_density = 1.81
drag_coefficent1 = .75
canopy_area1 = 15.029593683
max_height = 3500

terminal_velocity1 = math.sqrt((2*rocket_mass*gravity)/(canopy_area1*drag_coefficent1*air_density))
decent_time1 = max_height/terminal_velocity1

print('Terminal Velocity: ', terminal_velocity1, 'm/s')
print('Decent Time: ',decent_time1, 'seconds')
'''

recovery_bay_start = vehicle.rocket_dict_recovery["recovery_bay"]["bottom_distance_from_aft"]  # m
cg = loads.cg_max_q
inertia = loads.inertia

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

def calcLateralAcceleration(drag_force, gravity, total_mass):
    '''
    drag_force: Parachute drag force [N]
    gravity: Gravitational acceleration [m / s^2]
    total_mass: Total mass of rocket [kg]
    ay: Lateral acceleration [m / s^2]
    '''
    ay = (drag_force / total_mass) - gravity
    return ay

def calcAngularAcceleration(drag_force, recovery_bay_start, inertia, cg):
    '''
    drag_force: Parachute drag force [N]
    recovery_bay_start: Location of start of recovery bay [m]
    inertia: Rotational inertia around center of gravity [kg m^2]
    cg: Location of center of gravity of rocket [m]
    r: Angular acceleration [1 / s^2]
    '''
    r = ((-1) * drag_force * (abs(recovery_bay_start - cg))) / inertia
    return r

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
    #print(np.cumsum(linear_density_array * dx))
    # print(shear_array) # TEST
    #print(f"Lateral acceleration: {ay}")
    #print(f"Nose lift: {noseLift}") # TEST
    #print(f"Fin lift: {finLift}") # TEST
    shear_array[int(recovery_bay_start / dx) - 1:] += drag_force

    return shear_array

def calcBending(shear_array, length_along_rocket_linspace):
    '''
    shear_array: Array of shear forces across rocket length [N]
    length_along_rocket_linspace: Array of rocket lengths [m]
    bending_array: Array of bending forces across rocket length [N m]
    '''
    dy = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]
    bending_array = np.cumsum(shear_array) * dy
    return bending_array

terminal_velocity = calcTerminalVelocity(total_mass, gravity, drag_coefficent, air_density, canopy_area) # solving for velocity setting weight and Drag equal
drag_force = calcDragForce(drag_coefficent, air_density, 5, canopy_area) # Formula from NASA website
descent_time = max_height/terminal_velocity

print ('Terminal Velocity: ', terminal_velocity, 'm/s')
print ("Descent Time: ", descent_time, 'seconds')
print ('Drag Force: ', drag_force, 'N')

ay = calcLateralAcceleration(drag_force, gravity, total_mass) # Lateral acceleration
r = calcAngularAcceleration(drag_force, recovery_bay_start, inertia, cg) # Angular acceleration
shear_array = np.array(calcShear(drag_force, recovery_bay_start, ay, linear_density_array, length_along_rocket_linspace, r, cg)) # Shear force array
bending_array = np.array(calcBending(shear_array, length_along_rocket_linspace)) # Bending moment array
plt.plot(length_along_rocket_linspace, shear_array)
plt.title("Shear Force vs Length Along Rocket")
plt.xlabel("Length Along Rocket [m]")
plt.ylabel("Shear Force [N]")
plt.grid()
plt.show()