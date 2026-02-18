import os
from pathlib import Path
import scipy.io as sio 
import numpy as np
import builtins
from dataclasses import fields, asdict
import shear_bolted_joints
import subprocess

from vehicle_parameters import parameters, wet_mass_distribution
import vehicle_parameters_functions
import SFD.RDOF_2.rdof_v2
import print_filter

repository_root_path, _ = vehicle_parameters_functions.Get_Repository_Root_Path()
PSPL_ROCKET_A_file_path = repository_root_path / Path(f"vehicle_parameters.csv")
structural_loads_path = (repository_root_path / "Structures_Analysis" / "structural_loads.mat").resolve()
structures_analysis_path = (repository_root_path / "Structures_Analysis").resolve()


export_path_list = [PSPL_ROCKET_A_file_path]

try:
    Six_DoF_csv_file_path = (
        repository_root_path
        / ".."          # one directory up from PSPL_Rocket_A
        / "PSPL-6DOF"
        / "TheSixDoF"
        / "Inputs"
        / "Saved Rockets"
        / "FUCK_MATLAB"
        / "vehicle_parameters.csv"
    ).resolve()

    # write to CSV for 6DOF to read since 6DOF is in matlabese
    Six_DoF_csv_file_path.parent.mkdir(parents=True, exist_ok=True)
    export_path_list.append(Six_DoF_csv_file_path)
except:
    print("6DOF export failed")

vehicle_parameters_functions.ExportObjectToCSV(parameters, PSPL_ROCKET_A_file_path)
vehicle_parameters_functions.ExportObjectToCSV(parameters, Six_DoF_csv_file_path)
vehicle_parameters_functions.ExportObjectToCSV(wet_mass_distribution, "wet_mass_distribution")

vehicle_parameters_functions.ExportObjectToMat(wet_mass_distribution, structures_analysis_path / "wet_mass_distribution.mat")

print("")

structural_loads_script_file_path = (structures_analysis_path / "Full_Rocket_Analysis.m")

launched_by = os.getenv("LAUNCHED_BY")

if launched_by != "matlab":
    environment_dictionary = os.environ.copy()
    environment_dictionary["LAUNCHED_BY"] = "python"
    
    matlab_executable_path = r"C:\Program Files\MATLAB\R2025b\bin\matlab.exe"
    
    completed_process = subprocess.run(
                    ["matlab", "-batch", f"run('{structural_loads_script_file_path}')"], 
                    check=True, 
                    env=environment_dictionary,
                    capture_output=False,
                )

    print(completed_process.stdout)
    print(completed_process.stderr)
    completed_process.check_returncode()


# if not invoked_by_matlab:
# else:
#     print("Skipping MATLAB launch because main.py was invoked by MATLAB.")



structural_loads = vehicle_parameters_functions.load_matlab_struct_as_dataclass(structural_loads_path)
# print(f"structural_loads.upper_strut.max_compreession: {structural_loads.upper_strut.max_compression} N")

parameters.unfreeze()
parameters.upper_strut_max_load = max(structural_loads.upper_strut.max_compression, structural_loads.upper_strut.max_tension)
parameters.mid_strut_max_load = max(structural_loads.mid_strut.max_compression, structural_loads.mid_strut.max_tension)
parameters.lower_strut_max_load = max(structural_loads.lower_strut.max_compression, structural_loads.lower_strut.max_tension)

parameters.fuel_tank_max_load = max(structural_loads.fuel_tank.max_compression, structural_loads.fuel_tank.max_tension)
parameters.oxygen_tank_max_load = max(structural_loads.oxygen_tank.max_compression, structural_loads.oxygen_tank.max_tension)
parameters.copv_tube_max_load = max(structural_loads.copv_tube.max_compression, structural_loads.copv_tube.max_tension)
parameters.freeze()


if __name__ == "__main__":
    with print_filter.context_manager(print_everything=False, print_margins=True, print_titles=True):
        shear_bolted_joints.Calculate_Shear_Bolted_Joints()
        pass

