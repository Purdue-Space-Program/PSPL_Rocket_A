import os
from pathlib import Path
import scipy.io as sio 
import numpy as np

# sorry
import sys,pathlib,collections,importlib.abc; r=next(p for p in pathlib.Path(__file__).resolve().parents if p.name=="PSPL_Rocket_A");   m=collections.defaultdict(list); [m[f.stem.casefold()].append(f) for f in r.rglob("*.py") if f.stem.isidentifier() and f.name != "__init__.py" and not (set(f.relative_to(r).parts)&{".git",".venv","__pycache__","build","dist"})]; dup={k:v for k,v in m.items() if   len(v)>1}; sys.meta_path.insert(0,type("AmbiguousBareImportBlocker",(importlib.abc.MetaPathFinder,),{"find_spec":lambda   self,fullname,path=None,target=None: (_ for _ in ()).throw(ImportError(f"The import {fullname!r} could refer to any following packages: "+" ".join("\n\t" + str(p.relative_to(r)) for p in dup[fullname.casefold()])+f"\n\nSpecify which package it is by using the folder.\nFor example:\n\t'import SFD.{fullname}'")) if "." not in fullname and fullname.casefold() in dup else None})()); sys.path.insert(0,str(r)) if str(r) not   in sys.path else None; [sys.path.append(str(v[0].parent)) for k,v in m.items() if k not in dup and str(v[0].parent) not in sys.path]

import vehicle_parameters
import vehicle_parameters_functions
import print_filter
import constants as c

import six_DoF_caller
import rdof_v2
import structural_loads
import shear_bolted_joints


def vehicle_analysis():
    
    repository_root_path, _ = vehicle_parameters_functions.Get_Repository_Root_Path()

    PSPL_ROCKET_A_file_path = repository_root_path / Path(f"vehicle_parameters.csv")
    six_DoF_file_path = (repository_root_path / ".." / "PSPL-6DOF"/ "TheSixDoF").resolve()
    structures_analysis_file_path = (repository_root_path / "Structures_Analysis").resolve()
    
    parameters, wet_mass_distribution, dry_mass_distribution = vehicle_parameters.main()
    parameters = six_DoF_caller.main(parameters)
    parameters = rdof_v2.main(parameters)
    parameters, wet_mass_distribution = structural_loads.main(parameters, wet_mass_distribution)
    parameters = shear_bolted_joints.main(parameters)

    vehicle_parameters_functions.ExportObjectToCSV(parameters, PSPL_ROCKET_A_file_path)
    vehicle_parameters_functions.ExportObjectToCSV(wet_mass_distribution, "wet_mass_distribution")



def main():
    vehicle_analysis()

if __name__ == "__main__":
    with print_filter.context_manager(print_everything=True, print_margins=True, print_titles=True):
        main()
