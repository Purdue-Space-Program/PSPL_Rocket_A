from CoolProp.CoolProp import PropsSI
import numpy as np
import press_sim
import constants as c

def CvCrit(qdot, Gs, T, P1):
    '''
    qdot: mass flow rate (ft^3/min)
    Gs: Specific gravity of gas (relative to specific_gravity_of_air=1)
    T: Absolute Upstream Temperature (Rankine)
    P1: Inlet pressure (psia)
    cv: Flow coefficient
    '''
    cv = qdot / (0.471 * 22.67 * P1 * np.sqrt(1/(Gs * T))) # https://www.swagelok.com/downloads/webcatalogs/en/ms-06-84.pdf
    return cv

def KelvinToRankine(T):
    return T * 1.8

# COPV
copv_temp = press_sim.T_COPV * c.KELVIN2RANK # [R]
copv_pressure = press_sim.P_COPV * c.PA2PSI # [psia]
copv_density = PropsSI('D', 'T', copv_temp, 'P', copv_pressure, 'Nitrogen') # kg/m^3
copv_gs = copv_density / 1.225 # Specific gravity of gas (relative to specific_gravity_of_air=1)
copv_qdot = (press_sim.V_dot_ox_nom + press_sim.V_dot_fu_nom) * c.M32FT3 * 60 # [ft^3/min]
copv_cv = CvCrit(copv_qdot, copv_gs, copv_temp, copv_pressure)
print(f"COPV Temperature: {copv_temp:.2f} R")
print(f"COPV Pressure: {copv_pressure:.2f} psia")
print(f"COPV Density: {copv_density:.2f} kg/m^3")
print(f"COPV Specific Gravity: {copv_gs:.2f}")
print(f"COPV Mass Flow Rate: {copv_qdot:.2f} ft^3/min")
print(f"COPV Cv: {copv_cv:.2f}")