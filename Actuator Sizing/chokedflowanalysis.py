import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.optimize import fsolve
from CoolProp.CoolProp import PropsSI
from main import *
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *

##### PARAMETERS #####
piston_area = np.pi * (piston_diameter / 2)**2 # [m^2]
shaft_area = np.pi * (shaft_diameter / 2)**2 # [m^2]
cp = PropsSI('C', 'P', 1*ATM2PA, 'T', 293, 'air')
cv = PropsSI('CVMASS', 'P', 1*ATM2PA, 'T', 293, 'air')
gamma = cp/cv # ratio of specific heats (should be around 1.4 and doesnt vary too much with temperature)
R = cp-cv # gas constant (should be around 287 J/kg.K)
cd = 0.6 # sharp edged orifice
stroke_length = 0 # initializing a variable to track movement of the piston
time = 0
dt = 0.0001
v1 = 6 * IN2M * (piston_area - shaft_area) # assuming that the piston bore is 6 inches from the end before actuation)
p1 = 1*ATM2PA
p2 = 1*ATM2PA
t1 = 273 # [K]
t2 = 273 # [K]
force_history = []
time_history = []
time = 0

def calc_mdot(piston_area, rho, cd, pressure):
    # https://en.wikipedia.org/wiki/Discharge_coefficient
    dp = 1*ATM2PA - pressure
    mdot = cd * piston_area * np.sqrt(2 * rho * dp)
    return mdot

def calc_p2(p1, v1, v2, gamma):
    # https://www.grc.nasa.gov/www/k-12/airplane/compexp.html
    p2 = p1 * (v1 / v2)**gamma
    return p2

def calc_t2(t1, p1, p2, gamma):
    # https://www.grc.nasa.gov/www/k-12/airplane/compexp.html
    t2 = t1 * (p2 / p1)**(1 - 1 / gamma)
    return t2

def calc_counter_force(p2):
    counter_force = p2 * piston_area
    return counter_force

def mach_func(mach):
    return (piston_area * p1 / np.sqrt(t1)) * np.sqrt(gamma / R) * mach * (1 + ((gamma - 1) / 2) * mach**2)**(-1 * (gamma + 1) / (2 * (gamma - 1)))

def calc_dist_vol(F_net, rod_mass, piston_diameter, time, piston_velocity, dist_travelled, shaft_area):
    time_step = 0.0001
    dist_travelled = dist_travelled + piston_velocity * time_step + 0.5 * (F_net / rod_mass) * time_step**2
    piston_velocity_new = piston_velocity + (F_net * time_step) / (rod_mass)
    piston_velocity = piston_velocity_new
    volume_swept = dist_travelled * (np.pi * (piston_diameter / 2)**2 - shaft_area)
    time += time_step
    return dist_travelled, volume_swept, time, piston_velocity

def calc_F_net(piston_force, piston_seal_length, shaft_seal_length, piston_seal_area, shaft_seal_area, force_valve):
    fc_piston = (4 * piston_seal_length * M2IN) * LBF2N # assuming 4, worst case, for now
    fc_shaft = (4 * shaft_seal_length * M2IN) * LBF2N
    fh_piston = (18 * piston_seal_area * M22IN2) * LBF2N # from parker oring handbook figure 5-10, lowkey using it for u-cup :skull:
    fh_shaft = (18 * shaft_seal_area * M22IN2) * LBF2N
    friction_piston = fc_piston + fh_piston
    friction_shaft = fc_shaft + fh_shaft
    f_net = piston_force - force_valve - friction_piston * 2 - friction_shaft # two seals on piston, 1 on rod
    return f_net

i = 1
piston_velocity = 0
while stroke_length <= piston_stroke_length:
    #mach = fsolve(mach_func, 0) # 0 is an initial guess
    #rho = PropsSI('D', 'P', p2, 'T', t2, 'air')
    #mdot = calc_mdot(piston_area, rho, cd, p2)
    piston_force = (pressure - p1) * np.pi * ((piston_diameter**2) / 4)
    F_net = calc_F_net(piston_force, piston_seal_length, shaft_seal_length, piston_seal_area, shaft_seal_area, 0) # assuming no force by valve for now
    stroke_length, v2, time, piston_velocity = calc_dist_vol(F_net, rod_mass, piston_diameter, time, piston_velocity, stroke_length, shaft_area)
    p1 = p2
    p2 = calc_p2(p1, v1, v2, gamma)
    
    if t1 <= 0:
        raise ValueError("fuck")
    
    t1 = t2
    t2 = calc_t2(t1, p1, p2, gamma)
    v1 = v2
    i+=1
    force_history.append(piston_force)
    time_history.append(time)
plt.plot(force_history, time_history)
plt.xlabel("Force [N]")
plt.ylabel("Time")
plt.title("Force vs Time")
plt.show()