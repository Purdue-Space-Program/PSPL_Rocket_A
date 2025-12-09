import math
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import vehicle_parameters as vehicle
import sfd as sfd
import loads as loads

LBF2N = 4.44822  # Pounds force to Newtons
FT2M = 0.3048  # Feet to Meters

rocket_dict_dry = vehicle.rocket_dict_dry
cg = loads.cg_max_q
inertia = loads.inertia
parachute_mass = vehicle.parachute_mass  # [kg]
recovery_bay_start = rocket_dict_dry["recovery_bay"]["bottom_distance_from_aft"]  # [m]
max_q_velocity = vehicle.parameters.max_velocity  # [m / s]
AOA = loads.AOA  # [radians]
velocity = max_q_velocity * np.sin(AOA) # [m / s]


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
gravity = 9.81
air_density = 1.81
drag_coefficent = 2.2
canopy_area = (10 * FT2M / 2)**2 * np.pi # [m^2]
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

def calcLateralAcceleration(drag_force, total_mass):
    '''
    drag_force: Parachute drag force [N]
    total_mass: Total mass of rocket [kg]
    ay: Lateral acceleration [m / s^2]
    '''
    ay = (drag_force / total_mass)
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

terminal_velocity = calcTerminalVelocity(total_mass, gravity, drag_coefficent, air_density, canopy_area) # solving for velocity setting weight and Drag equal
drag_force = calcDragForce(drag_coefficent, air_density, velocity, canopy_area) # Formula from NASA website
descent_time = max_height/terminal_velocity

print ('Terminal Velocity: ', terminal_velocity, 'm/s')
print ("Descent Time: ", descent_time, 'seconds')
print ('Drag Force: ', drag_force, 'N')

ay = calcLateralAcceleration(drag_force, total_mass) # Lateral acceleration
r = calcAngularAcceleration(drag_force, recovery_bay_start, inertia, cg) # Angular acceleration
shear_array = np.array(calcShear(drag_force, recovery_bay_start, ay, linear_density_array, length_along_rocket_linspace, r, cg)) # Shear force array
bending_array = np.array(calcBending(shear_array, length_along_rocket_linspace)) # Bending moment array
axial_array = np.array(calcAxial(drag_force, linear_density_array, length_along_rocket_linspace)) # Axial forces array

# Plotting
N2LBS = 0.224809
M2FT = 3.28084

plt.figure()
plot_num = 1
for variable in ["shear_array", "bending_array", "axial_array"]:
    if variable == "shear_array":
        plot = shear_array * N2LBS
        ylabel = "Shear Force [lbs]"
        title = f"Shear Forces at Recovery"
    if variable == "bending_array":
        plot = bending_array * N2LBS
        ylabel = "Bending Moment [lbs-ft]"
        title = f"Bending Moments at Recovery"
    if variable == "axial_array":
        plot = axial_array * N2LBS
        ylabel = "Axial Force [lbs]"
        title = f"Axial Forces at Recovery"
    plt.subplot(1,3, plot_num)
    plt.plot(length_along_rocket_linspace * M2FT, plot)
    plt.title(title)
    plt.xlabel("Length from aft [ft]")
    plt.ylabel(ylabel)
    plt.grid()
    plot_num += 1
plt.show()



'''
plt.plot(length_along_rocket_linspace, bending_array)
plt.title("Bending Moment vs Length Along Rocket")
plt.xlabel("Length Along Rocket [m]")
plt.ylabel("Bending Moment [Nm]")
plt.grid()
plt.show()
'''