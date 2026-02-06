import os
from pathlib import Path
from scipy.io import savemat
from dataclasses import fields, asdict

import vehicle_parameters
import vehicle_parameters_functions
import SFD.rdof

######### paused until i figure ts out ##########
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
    pass


vehicle_parameters_functions.ConvertObjectToCSV(vehicle_parameters.parameters, "vehicle_parameters")
# vehicle_parameters_functions.ConvertObjectToCSV(vehicle_parameters.wet_mass_distribution, "wet_mass_distribution")


def convert_mass_distribution_to_matlab_dict(mass_distribution_object):
    matlab_struct_dictionary = {}
    for dataclass_field in fields(mass_distribution_object):
        mass_component_object = getattr(mass_distribution_object, dataclass_field.name)
        matlab_struct_dictionary[dataclass_field.name] = asdict(mass_component_object)
    return matlab_struct_dictionary

os.chdir("..")

savemat(
    r"Structures Analysis/wet_mass_distribution.mat",
    {"wet_mass_distribution": convert_mass_distribution_to_matlab_dict(vehicle_parameters.wet_mass_distribution)}
)