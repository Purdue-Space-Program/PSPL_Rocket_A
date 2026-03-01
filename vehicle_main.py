import os
from pathlib import Path
import scipy.io as sio 
import numpy as np
import subprocess

from vehicle_parameters import parameters, wet_mass_distribution
import vehicle_parameters_functions
import print_filter
import constants as c

import SFD.RDOF_2.rdof_v2 as rdof_v2
import Structures_Analysis.structural_loads as structural_loads
import shear_bolted_joints


def vehicle_analysis():
    repository_root_path, _ = vehicle_parameters_functions.Get_Repository_Root_Path()

    PSPL_ROCKET_A_file_path = repository_root_path / Path(f"vehicle_parameters.csv")
    six_DoF_file_path = (repository_root_path / ".." / "PSPL-6DOF"/ "TheSixDoF").resolve()
    structures_analysis_file_path = (repository_root_path / "Structures_Analysis").resolve()

    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    END = '\033[0m'
    GREEN = '\033[92m'

    vehicle_parameters_functions.ExportObjectToMat(wet_mass_distribution, structures_analysis_file_path / "wet_mass_distribution.mat")
    
    six_DoF_vehicle_parameters_csv_file_path = (
        six_DoF_file_path
        / "Inputs"
        / "Saved Rockets"
        / "FUCK_MATLAB"
        / "vehicle_parameters.csv"
    ).resolve()

    six_DoF_script_file_path = six_DoF_file_path / "non_GUI_run.m"

    # write to CSV for 6DOF to read since 6DOF is in matlabese
    vehicle_parameters_functions.ExportObjectToCSV(parameters, six_DoF_vehicle_parameters_csv_file_path)

    # weird thing where matlab was using venv version of python from Rocket_A    
    # fuck = subprocess.run(["33matlab", "-batch", "pyenv"])
        
    six_DoF_completed_process = None
    try:
        matlab_executable_path = r"matlab" #r"C:\Program Files\MATLAB\R2025b\bin\matlab.exe"
        # print(f"\n\n{matlab_executable_path}", always_print_this=True)

        startup_information = subprocess.STARTUPINFO()
        startup_information.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        startup_information.wShowWindow = subprocess.SW_HIDE

        creation_flags = subprocess.CREATE_NO_WINDOW | subprocess.DETACHED_PROCESS

        six_DoF_completed_process = subprocess.run(
            [matlab_executable_path, "-batch", f"run('{six_DoF_script_file_path}')"],
            check=True,
            capture_output=True,
            text=True,
            startupinfo=startup_information,
            creationflags=creation_flags,
        )

        six_DoF_output = vehicle_parameters_functions.load_matlab_struct_as_dataclass(
            six_DoF_file_path / "output.mat"
        )

    except subprocess.CalledProcessError as matlab_error:
        print(matlab_error.stdout, always_print_this=True)
        print(matlab_error.stderr, always_print_this=True)
        six_DoF_completed_process.check_returncode()
        raise
        print(f"{RED} =============6DOF run failed, using last saved outputs, take results with grain of salt============= {END}", always_print_this = True)

    except Exception as other_error:
        if six_DoF_completed_process is not None:
            print(six_DoF_completed_process.stdout, always_print_this=True)
            print(six_DoF_completed_process.stderr, always_print_this=True)
            six_DoF_completed_process.check_returncode()
        raise
        print(f"{RED} =============6DOF run failed, using last saved outputs, take results with grain of salt============= {END}", always_print_this = True)


    six_DoF_output = vehicle_parameters_functions.load_matlab_struct_as_dataclass(six_DoF_file_path / "output.mat")

    parameters.unfreeze()
    # parameters.six_DoF_off_the_rail_acceleration = six_DoF_output.
    # parameters.six_DoF_off_the_rail_velocity = six_DoF_output.
    # parameters.six_DoF_max_acceleration = six_DoF_output.
    parameters.six_DoF_max_Q_mach = six_DoF_output.max_Q_Mach
    parameters.six_DoF_max_Q_velocity = parameters.six_DoF_max_Q_mach * c.SPEED_OF_SOUND
    parameters.six_DoF_max_Q_horizontal_velocity = six_DoF_output.max_Q_horizontal_velocity 
    parameters.six_DoF_max_Q_acceleration = six_DoF_output.max_Q_acceleration
    parameters.six_DoF_apogee_horizontal_velocity = six_DoF_output.apogee_horizontal_velocity
    parameters.six_DoF_estimated_apogee = six_DoF_output.apogee_altitude
    print(f"parameters.six_DoF_estimated_apogee: {parameters.six_DoF_estimated_apogee:.2f} m")
    print(f"parameters.six_DoF_estimated_apogee: {parameters.six_DoF_estimated_apogee * c.M2FT:.2f} ft")
    # parameters.six_DoF_total_impulse = six_DoF_output.

    parameters.freeze()

    rdof_v2.main()
    structural_loads.main()
    shear_bolted_joints.main()

    vehicle_parameters_functions.ExportObjectToCSV(parameters, PSPL_ROCKET_A_file_path)
    vehicle_parameters_functions.ExportObjectToCSV(wet_mass_distribution, "wet_mass_distribution")


def main():
    vehicle_analysis()

if __name__ == "__main__":
    with print_filter.context_manager(print_everything=True, print_margins=True, print_titles=True):
        main()
