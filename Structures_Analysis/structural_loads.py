import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.ticker import MaxNLocator, MultipleLocator
from scipy.integrate import solve_ivp
from scipy.io import savemat
from pathlib import Path
import subprocess



def calculate_structural_loads():
    import sys
    import os

    os.chdir(os.path.dirname(__file__))
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

    from vehicle_parameters import parameters
    import vehicle_parameters_functions
    import vehicle_main
    import constants as c
    # print("\n\n")
    # print(f"os.getcwd: {os.getcwd()}")
    # print(f"sys.path: {sys.path}")
    # print("\n\n")
    
    repository_root_path, _ = vehicle_parameters_functions.Get_Repository_Root_Path()
    structures_analysis_file_path = (repository_root_path / "Structures_Analysis").resolve()
    structural_loads_file_path = (repository_root_path / "Structures_Analysis" / "structural_loads.mat").resolve()
    structural_loads_script_file_path = (structures_analysis_file_path / "Full_Rocket_Analysis.m")

    
    launched_by = os.getenv("LAUNCHED_BY")
    if launched_by != "matlab":
        environment_dictionary = os.environ.copy()
        environment_dictionary["LAUNCHED_BY"] = "python"
        
        matlab_executable_path = r"C:\Program Files\MATLAB\R2025b\bin\matlab.exe"
        
        structural_loads_completed_process = subprocess.run(
                        ["matlab", "-batch", f"run('{structural_loads_script_file_path}')"], 
                        check=True, 
                        env=environment_dictionary,
                        capture_output=False,
                    )

        print(structural_loads_completed_process.stdout)
        print(structural_loads_completed_process.stderr)
        structural_loads_completed_process.check_returncode()

    structural_loads = vehicle_parameters_functions.load_matlab_struct_as_dataclass(structural_loads_file_path)

    parameters.unfreeze()
    parameters.upper_strut_max_load = max(structural_loads.upper_strut.max_compression, structural_loads.upper_strut.max_tension)
    parameters.mid_strut_max_load = max(structural_loads.mid_strut.max_compression, structural_loads.mid_strut.max_tension)
    parameters.lower_strut_max_load = max(structural_loads.lower_strut.max_compression, structural_loads.lower_strut.max_tension)

    parameters.fuel_tank_max_load = max(structural_loads.fuel_tank.max_compression, structural_loads.fuel_tank.max_tension)
    parameters.oxygen_tank_max_load = max(structural_loads.oxygen_tank.max_compression, structural_loads.oxygen_tank.max_tension)
    parameters.copv_tube_max_load = max(structural_loads.copv_tube.max_compression, structural_loads.copv_tube.max_tension)
    parameters.freeze()

def main():
    calculate_structural_loads()

if __name__ == "__main__":
    main()