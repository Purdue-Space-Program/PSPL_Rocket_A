import vehicle_parameters
import vehicle_parameters_functions
# import SFD.loads
from pathlib import Path



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