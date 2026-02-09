starting_directory = pwd;
cd('C:..')

python_executable_path = "python";
python_script_path = "main.py";
[exit_status, command_output] = system(python_executable_path + " " + python_script_path);

% fprintf("exit_status: %d\n", exit_status);
% fprintf("command_output:\n%s\n", command_output);

vehicle_parameters_csv_path = fullfile("vehicle_parameters.csv");
csv_import_options = detectImportOptions(vehicle_parameters_csv_path);
csv_import_options.CommentStyle = "#";
csv_import_options.VariableNamesLine = 3;
csv_import_options.DataLines = [4, Inf];
vehicle_parameters_csv = readtable(vehicle_parameters_csv_path, csv_import_options);

% convert to struct for syntactic sugar
vehicle_parameters = cell2struct( ...
    num2cell(vehicle_parameters_csv.value), ...
    cellstr(vehicle_parameters_csv.parameter_name), ...
    1);

mat_file_path = fullfile("Structures Analysis\vehicle_parameters.mat");
save(mat_file_path, "vehicle_parameters", "vehicle_parameters_csv");

cd(starting_directory)