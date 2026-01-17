import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *
import vehicle_parameters as v

##### INPUT PARAMETERS #####
inner_diameter = v.parameters.tank_inner_diameter # [m]
wall_thickness = v.parameters.tank_wall_thickness # [m]
tank_pressure = v.parameters.tank_pressure # [Pa]
aluminum_6061_T6_yield_strength = 35000 * PSI2PA # [Pa]
aluminum_6061_T6_ultimate_strength = 42000 * PSI2PA # [Pa]
FoS_yield = 1.5
FoS_ultimate = 2


def calc_MoS(limit_load_stress, max_allowable_stress, FoS):
    design_load = limit_load_stress * FoS
    MoS = (max_allowable_stress / design_load) - 1
    return MoS

def calc_hoop_stress(pressure, inner_diameter, thickness):
    material_stress = pressure * (inner_diameter/2) / thickness
    return material_stress

limit_load_stress = calc_hoop_stress(tank_pressure, inner_diameter, wall_thickness)

MoS_yield = calc_MoS(limit_load_stress, aluminum_6061_T6_yield_strength, FoS_yield)
MoS_ultimate = calc_MoS(limit_load_stress, aluminum_6061_T6_ultimate_strength, FoS_ultimate)

print(f"\nYield FoS: {FoS_yield}, Ultimate FoS: {FoS_ultimate}")
print(f"Tank Pressure: {tank_pressure * PA2PSI:.2f} psi")

print(f"\nYield stress MoS: {MoS_yield:.3f}")
print(f"Ultimate stress MoS: {MoS_ultimate:.3f}")

