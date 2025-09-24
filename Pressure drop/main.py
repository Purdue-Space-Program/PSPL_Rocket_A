from fluids import fittings
from fluids import core
from CoolProp import CoolProp as CP
from pint import UnitRegistry
import math

u = UnitRegistry()

def find_drop_ox(K, type, line_velocity_ox):
    dP = core.dP_from_K(K, rho=rho_ox.magnitude, V=line_velocity_ox.magnitude) * u.Pa
    dP.to('psi')
    print(f'LOX LINE: Pressure Drop at {type}: {dP:.2f}')

def sharp_edged_inlet():
    K = 0.5
    return K

def inward_projecting_sharp_edged_inlet():
    K = 0.78
    return K

def sharp_edged_inlet_at_angle(angle):
    K = 0.5 + 0.3 * math.cos(angle) + 0.2 * math.cos(angle)**2
    return K

def rounded_inlet(radius, diameter):
    r_by_d = radius / diameter
    if r_by_d > 0.2:
        K = 0.03
    else:
        K = 0.5 - 1.07 * r_by_d**(0.5) - 2.13 * r_by_d + 8.24 * r_by_d**(1.5) - 8.48 * r_by_d**2 + 2.90 * r_by_d**(2.5)
    return K