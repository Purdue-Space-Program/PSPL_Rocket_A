from CoolProp.CoolProp import PropsSI
import numpy as np

def CvCrit(qdot, Gs, T, P1):
    '''
    qdot: mass flow rate (ft^3/min)
    Gs: Specific gravity of gas (relative to air=1)
    T: Absolute Upstream Temperature (Rankine)
    P1: Inlet pressure (psia)
    cv: Flow coefficient
    '''
    cv = qdot / (0.471 * 22.67 * P1 * np.sqrt(1/(Gs * T))) # https://www.swagelok.com/downloads/webcatalogs/en/ms-06-84.pdf
    return cv

def KelvinToRankine(T):
    return T * 1.8

