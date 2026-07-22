import math
from tabulate import tabulate

# sorry
import sys,pathlib,collections,importlib.abc; r=next(p for p in pathlib.Path(__file__).resolve().parents if p.name=="PSPL_Rocket_A");   m=collections.defaultdict(list); [m[f.stem.casefold()].append(f) for f in r.rglob("*.py") if f.stem.isidentifier() and f.name != "__init__.py" and not (set(f.relative_to(r).parts)&{".git",".venv","__pycache__","build","dist"})]; dup={k:v for k,v in m.items() if   len(v)>1}; sys.meta_path.insert(0,type("AmbiguousBareImportBlocker",(importlib.abc.MetaPathFinder,),{"find_spec":lambda   self,fullname,path=None,target=None: (_ for _ in ()).throw(ImportError(f"The import {fullname!r} could refer to any following packages: "+" ".join("\n\t" + str(p.relative_to(r)) for p in dup[fullname.casefold()])+f"\n\nSpecify which package it is by using the folder.\nFor example:\n\t'import SFD.{fullname}'")) if "." not in fullname and fullname.casefold() in dup else None})()); sys.path.insert(0,str(r)) if str(r) not   in sys.path else None; [sys.path.append(str(v[0].parent)) for k,v in m.items() if k not in dup and str(v[0].parent) not in sys.path]

import constants as c
import vehicle_parameters_functions 
import vehicle_parameters
import vehicle_main
import print_filter

def Calculate_Pressure(parameters): # this is a separate function for reasons i dont remember but its important...
    
    def Round_Up_To_Nearest_Multiple(original_value, multiple):
        rounded_value = multiple * math.ceil(original_value / multiple)
        return(rounded_value)
    
    
    parameters.unfreeze()
    parameters.maximum_tank_pressure_to_acount_for_droop = parameters.nominal_tank_pressure / 0.75 # the set pressure is higher than the nominal tank pressure to account for droop that occurs when the run valves open and there is flow through the regulator
    
    relief_valve_set_pressure = Round_Up_To_Nearest_Multiple(original_value = parameters.maximum_tank_pressure_to_acount_for_droop, multiple = (25 * c.PSI2PA)) # relief valves sold in increments of 25 psi
    
    parameters.maximum_allowable_tank_pressure = relief_valve_set_pressure * 1.1 # relief valves are rated for full flow when they are pressurized at 10% over their set pressure
    parameters.hydroproof_tank_pressure = parameters.maximum_allowable_tank_pressure * parameters.hydroproof_factor # this is the exact pressure to hydroproof at to prove that the tanks can handle the pressure that the relief valves will allow the tank to be
    parameters.freeze()

    output_table_data = []
    
    output_table_data.append(("nominal_tank_pressure", f"{parameters.nominal_tank_pressure * c.PA2PSI:.2f} PSI"))
    output_table_data.append(("maximum_tank_pressure_to_acount_for_droop", f"{parameters.maximum_tank_pressure_to_acount_for_droop * c.PA2PSI:.2f} PSI"))
    output_table_data.append(("relief_valve_set_pressure", f"{relief_valve_set_pressure * c.PA2PSI:.2f} PSI"))
    output_table_data.append(("maximum_allowable_tank_pressure", f"{parameters.maximum_allowable_tank_pressure * c.PA2PSI:.2f} PSI"))
    output_table_data.append(("hydroproof_tank_pressure", f"{parameters.hydroproof_tank_pressure * c.PA2PSI:.2f} PSI"))
    
    print(tabulate(output_table_data, tablefmt = "simple_grid"))
    
    return(parameters)
    
    
    

def main(parameters):
    parameters = Calculate_Pressure(parameters)
    return(parameters)

if __name__ == "__main__":
    parameters, wet_mass_distribution, dry_mass_distribution = vehicle_parameters.main()
    main(parameters)
