import os
from pathlib import Path
import scipy.io as sio 
import numpy as np
import builtins
from dataclasses import fields, asdict
import shear_bolted_joints

from vehicle_parameters import parameters, wet_mass_distribution
import vehicle_parameters_functions
# import SFD.rdof
import print_filter

repository_root_path, _ = vehicle_parameters_functions.Get_Repository_Root_Path()

PSPL_ROCKET_A_file_path = repository_root_path / Path(f"vehicle_parameters.csv")

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




# sio.savemat(
#     structures_analysis_path / "wet_mass_distribution.mat",
#     {"wet_mass_distribution": convert_mass_distribution_to_matlab_dict(vehicle_parameters.wet_mass_distribution)}
# )

# invoked_by_matlab = os.getenv("PSPL_INVOKED_BY_MATLAB") == "1"

#     ], check=True)
# if not invoked_by_matlab:
#     import subprocess
#     full_rocket_analysis_script = (structures_analysis_path / "Full_Rocket_Analysis.m").as_posix()
#     subprocess.run([
#         "matlab",
#         "-batch",
#         f"run('{full_rocket_analysis_script}')"
#     ], check=True)
# else:
#     print("Skipping MATLAB launch because main.py was invoked by MATLAB.")





structural_loads_path = (repository_root_path / "Structures_Analysis" / "structural_loads.mat").resolve()
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

print(f"parameters.upper_strut_max_load: {parameters.upper_strut_max_load:.2f}")


if __name__ == "__main__":
    with print_filter.context_manager(print_everything=False, print_margins=True, print_titles=True):
        shear_bolted_joints.Calculate_Shear_Bolted_Joints()

