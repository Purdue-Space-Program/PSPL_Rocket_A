import csv
from pathlib import Path
import sys
import subprocess
from datetime import datetime
import inspect
import io
from dataclasses import dataclass, fields, field, make_dataclass, asdict
import scipy.io as sio 
import numpy as np
import vehicle_parameters


def load_matlab_struct_as_dataclass(file_path_string):
    matlab_struct_name = Path(file_path_string.parts[-1]).stem
    weird_matlab_struct = sio.loadmat(file_path_string, struct_as_record=False, squeeze_me=True)
    normal_matlab_data_struct = weird_matlab_struct[matlab_struct_name]
    
    # for field_name in normal_matlab_data_struct._fieldnames:
    #     print(f"field_name: {field_name}")    
    
    
    # the_dataclass = make_dataclass(
    #     matlab_struct_name,
    #     [(normal_matlab_data_struct._fieldnames, type(field_values_dictionary[field_name])) for field_name in field_values_dictionary]
    # )
    # @dataclass(frozen=True)
    # class MassComponent:
    # dataclass = {}


        # clean_dictionary[key] = Convert_Matlab_Struct_To_Python_Dictionary(
        #     raw_matlab_data_dictionary[key]
        # )

    return normal_matlab_data_struct


def Convert_Matlab_Struct_To_Python_Dictionary(matlab_object):
    if hasattr(matlab_object, "_fieldnames"):
        field_values_dictionary = {}

        for field_name in matlab_object._fieldnames:
            field_values_dictionary[field_name] = Convert_Matlab_Struct_To_Python_Dictionary(
                getattr(matlab_object, field_name)
            )

        dataclass_type = make_dataclass(
            "MatlabStruct",
            [(field_name, type(field_values_dictionary[field_name])) for field_name in field_values_dictionary]
        )

        return dataclass_type(**field_values_dictionary)


def _convert_struct_array(matlab_struct_array):
    python_dictionary = {}

    struct_element = matlab_struct_array[0, 0]

    for field_name in struct_element.dtype.names:
        field_value = struct_element[field_name]
        python_dictionary[field_name] = Convert_Matlab_Struct_To_Python_Dictionary(field_value)

    return python_dictionary




def convert_mass_distribution_to_matlab_dict(mass_distribution_object):
    matlab_struct_dictionary = {}
    for dataclass_field in fields(mass_distribution_object):
        mass_component_object = getattr(mass_distribution_object, dataclass_field.name)
        matlab_struct_dictionary[dataclass_field.name] = asdict(mass_component_object)
    return matlab_struct_dictionary



def ExportObjectToCSV(object, export_file_path):
    repository_root_path, caller_file_path = Get_Repository_Root_Path()
    export_file_path = Path(export_file_path)

    if export_file_path.suffix != ".csv":
        export_file_path = export_file_path.with_suffix(".csv")
    # if not export_file_path.is_absolute():
    #     export_file_path = repository_root_path / export_file_path
    

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
        
        csv_writer_handle.writerow(["parameter_name", "value"])
        
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

            print(f"Mass Distribution CSV Exported to {export_file_path}")

        elif isinstance(object, vehicle_parameters.VehicleParameters):
            vehicle_parameters_object = object
            
            for field_object in fields(vehicle_parameters_object):
                if field_object.name.startswith("_"):
                    continue
                csv_writer_handle.writerow([field_object.name, getattr(vehicle_parameters_object, field_object.name)])
            
            print(f"Vehicle Parameters CSV Exported to {export_file_path}")
    
        else:
            raise ValueError("da fuq")








def ExportObjectToMat(object, export_file_path):
    sio.savemat(
        export_file_path,
        {"wet_mass_distribution": convert_mass_distribution_to_matlab_dict(object)}
    )





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
