import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from CoolProp.CoolProp import PropsSI
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *
import vehicle_parameters as v

##### OPERATING CONDITIONS #####
g = GRAVITY
n2_temperature = 0 # [K] Temperature of N2 at diffuser exit???
n2_mdot = 0 # [kg/s]
ullage_pressure = v.parameters.tank_pressure # [Pa] Diffuser exit pressure of N2 / Ullage Pressure
temperature_wall = 0
dTemp = abs(temperature_wall-n2_temperature)
tank_height = v.parameters.fuel_tank_length # [m]
tank_diameter = v.parameters.tank_inner_diameter # [m]

##### INPUT PARAMETERS #####
desired_velocity = 1 # [m/s]
cd = 0.65 # Conservative estimate

def calc_diffuser_area(n2_temperature, ullage_pressure, n2_mdot, desired_velocity, cd):
    rho = PropsSI('D', 'P', ullage_pressure, 'T', n2_temperature, 'Nitrogen')
    diffuser_area = n2_mdot / (rho * desired_velocity * cd)
    return diffuser_area

def calc_ri(n2_temperature, g, charLength, dTemp, desired_velocity):
    cte = 1 / n2_temperature
    ri = cte * g * charLength * dTemp / desired_velocity**2
    return ri

def check_choked(desired_velocity, n2_temperature, ullage_pressure):
    speed_of_sound = PropsSI('A', 'P', ullage_pressure, 'T', n2_temperature, 'Nitrogen')
    M = desired_velocity / speed_of_sound
    return M

if __name__ == "__main__":
    
    if tank_height >= tank_diameter:
    

    charLength = 0

    M = check_choked(desired_velocity, n2_temperature, ullage_pressure)
    if M >=1:
        print(f"Mach Number: {M}. Flow is choked!!!!! Ullage Collapse!!!!!!!!")
    elif M > 0.3:
        print(f"Mach Number: {M}. Outside of safe range!!!!! Reconsider Parameters!!!!!")
    elif M >= 0:
        print(f"Mach Number: {M}. Flow is (probably) not choked!!!!!!")
    else:
        print(f"Mach Number: {M}. Invalid Output!!!!! Check Parameters!!!!!")