import numpy as np
from CoolProp.CoolProp import PropsSI

from old_CvCalc import Gs
def mass_flow_rate(Cv, P, T, N2=22.57, Gs=0.967):
    '''
    Cv: Flow coefficient
    P: Inlet pressure (psia)
    T: Absolute Upstream Temperature (Rankine)
    N2: Constant for mass flow units (default is 22.57 for scfm)
    Gs: Specific gravity of gas (relative to specific_gravity_of_air=1)
    qdot: mass flow rate (ft^3/min)
    '''
    qdot = 0.471 * N2 * Cv * P * np.sqrt(1/(Gs * T)) # https://www.swagelok.com/downloads/webcatalogs/en/ms-06-84.pdf (page 3)
    return qdot

def CvCrit(qdot, P, T, N2=22.57, Gs=0.967):
    '''
    qdot: mass flow rate (ft^3/min)
    P: Inlet pressure (psia)
    T: Absolute Upstream Temperature (Rankine)
    N2: Constant for mass flow units (default is 22.57 for scfm)
    Gs: Specific gravity of gas (relative to specific_gravity_of_air=1)
    cv: Flow coefficient
    '''
    cv = qdot / (0.471 * N2 * P * np.sqrt(1/(Gs * T))) # https://www.swagelok.com/downloads/webcatalogs/en/ms-06-84.pdf (page 3)
    return cv

rho_n2 = PropsSI('D', 'T', 293.15, 'P', 101325, 'Nitrogen') # kg/m^3

# Regulator
Cv = 0.8
P = 3300 # psia
T = 293.15 * 1.8 + 32 # R

qdot = mass_flow_rate(Cv, P, T)
print(f"Mass Flow Rate: {qdot:.2f} ft^3/min")