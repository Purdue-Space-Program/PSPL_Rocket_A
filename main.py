import os
from pathlib import Path
import scipy.io as sio 
import numpy as np
import builtins
from dataclasses import fields, asdict
import shear_bolted_joints

import vehicle_parameters
import vehicle_parameters_functions
import SFD.rdof
import print_filter


def Convert_Matlab_Struct_To_Python_Dictionary(matlab_object):
    """
    Recursively convert scipy.io.loadmat struct objects
    into native Python dictionaries.
    """

    if isinstance(matlab_object, np.ndarray):

        # If it is a MATLAB struct array
        if matlab_object.dtype.names is not None:
            return _convert_struct_array(matlab_object)

        # If it is a numeric array
        if matlab_object.size == 1:
            return matlab_object.item()

        return matlab_object

    return matlab_object


def _convert_struct_array(matlab_struct_array):
    python_dictionary = {}

    struct_element = matlab_struct_array[0, 0]

    for field_name in struct_element.dtype.names:
        field_value = struct_element[field_name]
        python_dictionary[field_name] = Convert_Matlab_Struct_To_Python_Dictionary(field_value)

    return python_dictionary


def load_matlab_file_as_clean_dictionary(file_path_string):
    raw_matlab_data_dictionary = sio.loadmat(file_path_string, struct_as_record=False, squeeze_me=False)

    clean_dictionary = {}

    for key in raw_matlab_data_dictionary:
        if key.startswith("__"):
            continue

        clean_dictionary[key] = Convert_Matlab_Struct_To_Python_Dictionary(
            raw_matlab_data_dictionary[key]
        )

    return clean_dictionary




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
    print("6DOF export failed")
    pass


vehicle_parameters_functions.ConvertObjectToCSV(vehicle_parameters.parameters, "vehicle_parameters")
# vehicle_parameters_functions.ConvertObjectToCSV(vehicle_parameters.wet_mass_distribution, "wet_mass_distribution")


def convert_mass_distribution_to_matlab_dict(mass_distribution_object):
    matlab_struct_dictionary = {}
    for dataclass_field in fields(mass_distribution_object):
        mass_component_object = getattr(mass_distribution_object, dataclass_field.name)
        matlab_struct_dictionary[dataclass_field.name] = asdict(mass_component_object)
    return matlab_struct_dictionary

structures_analysis_path = (repository_root_path.parent / "Structures_Analysis").resolve()

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




# upper_strut_struct = structural_loads["upper_strut"][0][0]

# upper_strut_max_compression = upper_strut_struct["max_compression"][0][0]
# upper_strut_max_tension = upper_strut_struct["max_tension"][0][0]

# print("Upper Strut Max Compression:", upper_strut_max_compression)
# print("Upper Strut Max Tension:", upper_strut_max_tension)






















def Convert_Matlab_Struct_To_Python_Dictionary(matlab_object):
    """
    Recursively convert scipy.io.loadmat struct objects
    into native Python dictionaries.
    """

    if isinstance(matlab_object, np.ndarray):

        # If it is a MATLAB struct array
        if matlab_object.dtype.names is not None:
            return _convert_struct_array(matlab_object)

        # If it is a numeric array
        if matlab_object.size == 1:
            return matlab_object.item()

        return matlab_object

    return matlab_object


def _convert_struct_array(matlab_struct_array):
    python_dictionary = {}

    struct_element = matlab_struct_array[0, 0]

    for field_name in struct_element.dtype.names:
        field_value = struct_element[field_name]
        python_dictionary[field_name] = Convert_Matlab_Struct_To_Python_Dictionary(field_value)

    return python_dictionary


def load_matlab_file_as_clean_dictionary(file_path_string):
    raw_matlab_data_dictionary = sio.loadmat(file_path_string, struct_as_record=False, squeeze_me=False)

    clean_dictionary = {}

    for key in raw_matlab_data_dictionary:
        if key.startswith("__"):
            continue

        clean_dictionary[key] = Convert_Matlab_Struct_To_Python_Dictionary(
            raw_matlab_data_dictionary[key]
        )

    return clean_dictionary


if not invoked_by_matlab:
    structural_loads = Convert_Matlab_Struct_To_Python_Dictionary(sio.loadmat(structures_analysis_path / "structural_loads.mat")["structural_loads"])
else:
    structural_loads = None

# print(structural_loads)

if __name__ == "__main__":
    with print_filter.context_manager(print_everything=False, print_margins=True, print_titles=True):
        shear_bolted_joints.Calculate_Shear_Bolted_Joints()
