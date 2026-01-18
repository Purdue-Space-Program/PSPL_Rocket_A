import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from CoolProp.CoolProp import PropsSI
from main import *
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *

##### PARAMETERS #####
piston_area = np.pi * (piston_diameter / 2)**2 # [m^2]
cp = PropsSI('C', 'P', 1*ATM2PA, 'T', 293, 'air')
cv = PropsSI('CVMASS', 'P', 1*ATM2PA, 'T', 293, 'air')
gamma = cp/cv # ratio of specific heats (should be around 1.4 and doesnt vary too much with temperature)
print(gamma)
R = cp-cv # gas constant (should be around 287 J/kg.K)
cd = 0.6

def calc_mdot(piston_area, rho, dp, pressure):
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

