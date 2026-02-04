import csv
import vehicle_parameters
from pathlib import Path
import sys
import subprocess
from datetime import datetime
import inspect


def ConvertObjectToCSV(object, file_name):
    
    export_file_path = Path(f"{file_name}.csv")
    
    try:
        caller_file_path = Path(selected_path.relative_to(repository_root_path).as_posix())
    except Exception:
        caller_file_path = Path(selected_path.as_posix())
    # print(f"caller path: {caller_file_path}")

    timestamp_string = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    
    with open(export_file_path, "w", newline="", encoding="utf-8") as csv_file_handle:
        
        csv_writer_handle = csv.writer(csv_file_handle)

        # fuck epoch
        csv_file_handle.write(f"# Accessed: {timestamp_string} (format: YYYY-MM-DD_HH-MM-SS)\n")
        csv_file_handle.write(f"# Accessed by: {caller_file_path.as_posix()}\n")
        
        

        if isinstance(object, vehicle_parameters.MassDistribution):
            
            csv_writer_handle.writerow([
                "name",
                "mass",
                "bottom_distance_from_aft",
                "length",
                "top_distance_from_aft",
            ])

            for component in object:
                csv_writer_handle.writerow([
                    component.name,
                    component.mass,
                    component.bottom_distance_from_aft,
                    component.length,
                    component.top_distance_from_aft,
                ])



def export_to_csv(export_file_path):
    

    
    
    with open(export_file_path, "w", newline="") as csv_file_handle:
        

        csv_writer_handle = csv.writer(csv_file_handle)
        csv_writer_handle.writerow(["parameter_name", "value"])

        for field_object in fields(parameters):
            if field_object.name.startswith("_"):
                continue
            csv_writer_handle.writerow([field_object.name, getattr(parameters, field_object.name)])
        
        # print(f"Vehicle Parameters CSV Exported to {export_file_path}")
# print("")



def Generate_CSV_Bytes_From_Parameters(parameters):
    string_buffer = io.StringIO()
    csv_writer_handle = csv.writer(string_buffer)

    csv_writer_handle.writerow(["parameter_name", "value"])

    for field_object in fields(parameters):
        if field_object.name.startswith("_"):
            continue
        csv_writer_handle.writerow([
            field_object.name,
            getattr(parameters, field_object.name),
        ])

    return string_buffer.getvalue().encode("utf-8")


def read_csv_bytes_without_comments(csv_path_file_path, comment_prefix="#"):
    filtered_lines = []

    with open(csv_path_file_path, "rb") as file_handle:
        for raw_line in file_handle:
            stripped_line = raw_line.lstrip()
            if not stripped_line.startswith(comment_prefix.encode()):
                filtered_lines.append(raw_line)

    return b"".join(filtered_lines)

def Determine_if_CSV_Files_are_Equal(csv_file_path_a, csv_file_path_b):
    bytes_a = read_csv_bytes_without_comments(csv_file_path_a)
    bytes_b = read_csv_bytes_without_comments(csv_file_path_b)
    CSV_files_are_equal = bytes_a == bytes_b
    return CSV_files_are_equal







def

    # for recording what file accesed this script
    main_module = sys.modules.get("__main__")
    if getattr(main_module, "__file__", None):
        selected_path = Path(main_module.__file__).resolve()
    else:
        selected_path = None
        for frame_info in inspect.stack()[1:]:
            candidate_filename = frame_info.filename
            if candidate_filename and not candidate_filename.startswith("<"):
                candidate_path = Path(candidate_filename)
                if candidate_path.suffix == ".py":
                    selected_path = candidate_path.resolve()
                    break
        if selected_path is None:
            selected_path = Path.cwd().resolve()

    repository_root_path = Path(
        subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            cwd=selected_path.parent,
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
    ).resolve()

    python_file_dir = Path(__file__).resolve().parent

    # also output here for record keeping
    PSPL_ROCKET_A_records_file_path = repository_root_path / Path("vehicle_parameters_records")
    PSPL_ROCKET_A_records_file_path.mkdir(exist_ok=True)



    PSPL_ROCKET_A_new_record_file_path = PSPL_ROCKET_A_records_file_path / f"vehicle_parameters_{timestamp_string}.csv"

    PSPL_ROCKET_A_file_path = repository_root_path / f"vehicle_parameters.csv"

    export_path_list = [PSPL_ROCKET_A_file_path, PSPL_ROCKET_A_new_record_file_path]

    try:
        Six_DoF_csv_file_path = (
            python_file_dir
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
