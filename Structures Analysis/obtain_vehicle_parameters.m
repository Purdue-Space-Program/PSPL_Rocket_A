currentDirectory = pwd;
cd('C:..')

python_executable_path = "python";
python_script_path = "vehicle_parameters.py";
[exit_status, command_output] = system(python_executable_path + " " + python_script_path);

%disp(exit_status)
%disp(command_output)

vehicle_parameters_csv_path = fullfile("vehicle_parameters.csv");
vehicle_parameters_csv = readtable(vehicle_parameters_csv_path);

mat_file_path = fullfile("C:Structures Analysis\vehicle_parameters.mat");

save(mat_file_path, "vehicle_parameters_csv")

cd(currentDirectory)

% convert to struct for syntactic sugar
vehicle_parameters = cell2struct( ...
    num2cell(vehicle_parameters_csv.value), ...
    cellstr(vehicle_parameters_csv.parameter_name), ...
    1);

