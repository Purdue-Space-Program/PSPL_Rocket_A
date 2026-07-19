repository_directory_name = "PSPL_Rocket_A";
starting_directory = pwd;

[~, current_directory_name] = fileparts(pwd);

fuck = 0;
while ~strcmp(current_directory_name, repository_directory_name)
    fuck = fuck + 1;
    if fuck > 20
        cd(starting_directory);
        error("Couldn't Find Repository directory starting_directory: %s", starting_directory);
    end

    cd('..');
    [~, current_directory_name] = fileparts(pwd);
end


python_executable_file_path = "python";
python_script_file_path = "Vehicle_Level/vehicle_main.py";
vehicle_parameters_csv_file_path = fullfile("Vehicle_Level/vehicle_parameters.csv");
vehicle_parameters_mat_file_path = fullfile("Structures/Structural_Loads/vehicle_parameters.mat");


launched_by = getenv("LAUNCHED_BY");
% fprintf("launched_by: %s\n", launched_by);

if ~(strcmp(launched_by, "python"))
    setenv("LAUNCHED_BY", "matlab");
    [exit_status, command_output] = system(python_executable_file_path + " " + python_script_file_path);
    fprintf("exit_status: %d\n", exit_status);
    fprintf("command_output:\n%s\n", command_output);
end


csv_import_options = detectImportOptions(vehicle_parameters_csv_file_path);
csv_import_options.CommentStyle = "#";
csv_import_options.VariableNamesLine = 3;
csv_import_options.DataLines = [4, Inf];
vehicle_parameters_csv = readtable(vehicle_parameters_csv_file_path, csv_import_options);

% convert to struct for syntactic sugar
vehicle_parameters = cell2struct( ...
    num2cell(vehicle_parameters_csv.value), ...
    cellstr(vehicle_parameters_csv.parameter_name), ...
    1);

save(vehicle_parameters_mat_file_path, "vehicle_parameters", "vehicle_parameters_csv");

cd(starting_directory)
