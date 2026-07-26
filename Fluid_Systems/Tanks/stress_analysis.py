import numpy as np
import matplotlib.pyplot as plt
import sys
import sys,pathlib,collections,importlib.abc; r=next(p for p in pathlib.Path(__file__).resolve().parents if p.name=="PSPL_Rocket_A");   m=collections.defaultdict(list); [m[f.stem.casefold()].append(f) for f in r.rglob("*.py") if f.stem.isidentifier() and f.name != "__init__.py" and not (set(f.relative_to(r).parts)&{".git",".venv","__pycache__","build","dist"})]; dup={k:v for k,v in m.items() if   len(v)>1}; sys.meta_path.insert(0,type("AmbiguousBareImportBlocker",(importlib.abc.MetaPathFinder,),{"find_spec":lambda   self,fullname,path=None,target=None: (_ for _ in ()).throw(ImportError(f"The import {fullname!r} could refer to any following packages: "+" ".join("\n\t" + str(p.relative_to(r)) for p in dup[fullname.casefold()])+f"\n\nSpecify which package it is by using the folder.\nFor example:\n\t'import SFD.{fullname}'")) if "." not in fullname and fullname.casefold() in dup else None})()); sys.path.insert(0,str(r)) if str(r) not   in sys.path else None; [sys.path.append(str(v[0].parent)) for k,v in m.items() if k not in dup and str(v[0].parent) not in sys.path]

from constants import *
import Vehicle_Level.vehicle_parameters as v

##### INPUT PARAMETERS #####
inner_diameter = v.parameters.tank_inner_diameter # [m]
wall_thickness = v.parameters.tank_wall_thickness # [m]
tank_pressure = 930 * PSI2PA # [Pa]
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

print(f"Tank Hoop Stress: {limit_load_stress * PA2PSI / 1000:.1f} ksi")
print(f"\nYield stress MoS: {MoS_yield:.3f}")
print(f"Ultimate stress MoS: {MoS_ultimate:.3f}")
print(f"Zero margin yield pressure: {aluminum_6061_T6_yield_strength / FoS_yield * (2 * wall_thickness / inner_diameter) * PA2PSI:.1f} psi")
print(f"Tank explodes at {aluminum_6061_T6_ultimate_strength / FoS_ultimate * (2 * wall_thickness / inner_diameter) * PA2PSI:.1f} psi")