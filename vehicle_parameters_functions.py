import csv
from pathlib import Path
import sys
import subprocess
from datetime import datetime
import inspect
import io
from dataclasses import dataclass, fields, field

import vehicle_parameters

def ConvertObjectToCSV(object, file_name):
    repository_root_path, caller_file_path = Get_Repository_Root_Path()
    
    export_file_path = Path(f"{file_name}.csv")
    
    try:
        caller_file_path = Path(caller_file_path.relative_to(repository_root_path).as_posix())
    except Exception:
        caller_file_path = Path(caller_file_path.as_posix())
    # print(f"caller path: {caller_file_path}")

    timestamp_string = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    
    with open(export_file_path, "w", newline="", encoding="utf-8") as csv_file_handle:
        
        csv_writer_handle = csv.writer(csv_file_handle)

        # fuck epoch
        csv_file_handle.write(f"# Accessed: {timestamp_string} (format: YYYY-MM-DD_HH-MM-SS)\n")
        csv_file_handle.write(f"# Accessed by: {caller_file_path.as_posix()}\n")
        
        if isinstance(object, vehicle_parameters.MassDistribution):
            mass_distribution_object = object
            
            # csv_writer_handle.writerow([
            #     "name",
            #     "mass",
            #     "bottom_distance_from_aft",
            #     "length",
            #     "top_distance_from_aft",
            # ])

            # for component in mass_distribution_object:
            #     csv_writer_handle.writerow([
            #         component.name,
            #         component.mass,
            #         component.bottom_distance_from_aft,
            #         component.length,
            #         component.top_distance_from_aft,
            #     ])
            
            for field_object in fields(mass_distribution_object):
                if field_object.name.startswith("_"):
                    continue
                csv_writer_handle.writerow([field_object.name, getattr(mass_distribution_object, field_object.name)])


        elif isinstance(object, vehicle_parameters.VehicleParameters):
            vehicle_parameters_object = object
            
            for field_object in fields(vehicle_parameters_object):
                if field_object.name.startswith("_"):
                    continue
                csv_writer_handle.writerow([field_object.name, getattr(vehicle_parameters_object, field_object.name)])
            
            print(f"Vehicle Parameters CSV Exported to {export_file_path}")
    
        else:
            raise ValueError("da fuq")         
        print("")



# def Generate_CSV_Bytes_From_Parameters(parameters):
#     string_buffer = io.StringIO()
#     csv_writer_handle = csv.writer(string_buffer)

#     csv_writer_handle.writerow(["parameter_name", "value"])

#     for field_object in fields(parameters):
#         if field_object.name.startswith("_"):
#             continue
#         csv_writer_handle.writerow([
#             field_object.name,
#             getattr(parameters, field_object.name),
#         ])

#     return string_buffer.getvalue().encode("utf-8")






def Generate_CSV_Bytes_From_Class_Object(object):
    string_buffer = io.StringIO()
    csv_writer_handle = csv.writer(string_buffer)

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
    else:
        raise TypeError("Unsupported class object type for CSV byte generation.")

    return string_buffer.getvalue().encode("utf-8")


# def Generate_CSV_Bytes_From_Object_CSV(csv_path_file_path, comment_prefix="#"):
#     return read_csv_bytes_without_comments(
#         csv_path_file_path,
#         comment_prefix=comment_prefix,
#     )





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



def Get_Repository_Root_Path():
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
    
    return(repository_root_path, selected_path)
