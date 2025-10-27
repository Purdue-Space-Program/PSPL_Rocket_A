from fluids import fittings
from fluids import core
from CoolProp import CoolProp as CP
from pint import UnitRegistry
import math
import constants

u = UnitRegistry()

# Define fluid properties
fuel_type = "Isopropanol"
ox_type = "Oxygen"

tank_pressure_ox = 250 * u.psi
saturation_temp_ox = CP.PropsSI('T', 'P', tank_pressure_ox.magnitude, 'Q', 0, 'oxygen') * u.kelvin

rho_fuel = constants.DENSITY_IPA * (u.kilogram / u.meter**3)
rho_ox = CP.PropsSI("D", "P", tank_pressure_ox.magnitude + 10, "T", saturation_temp_ox.magnitude, ox_type) * (u.kilogram / u.meter**3)

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

def sharp_edged_sudden_contraction(diameter_one, diameter_two):
    K = 0.5 * (1 - (diameter_one / diameter_two)**2)
    return K

def sharp_edged_sudden_expansion(diameter_one, diameter_two):
    K = 1.05 * (1 - (diameter_one / diameter_two)**2)**2
    return K

def outlet():
    K = 1.05
    return K

def bend(angle, radius, diameter, friction_factor):
    r_by_d = radius / diameter
    K = friction_factor * angle * r_by_d + (0.1 + 2.4 * friction_factor) * math.sin(angle / 2) + (6.6 * friction_factor * math.sin(angle / 2)**0.5 + math.sin(angle / 2)) / (r_by_d**((4 * angle) / math.pi))
    return K

def valve(flow_coefficient, diameter):
    K = (890.4 * diameter**4) / (flow_coefficient**2)
    return K

find_drop_ox(sharp_edged_inlet(), "Sharp Edged Inlet", )