clear;
clc;

structural_loads = struct();
make_components


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

repository_path = pwd;
structural_loads_path = fullfile(repository_path, "Structures", "Structural_Loads");
cd(structural_loads_path);


%% Struts

[structural_loads.upper_strut.max_compression, structural_loads.upper_strut.max_tension] = StructuresAnalysis(upper_strut);
[structural_loads.mid_strut.max_compression, structural_loads.mid_strut.max_tension] = StructuresAnalysis(mid_strut);
[structural_loads.lower_strut.max_compression, structural_loads.lower_strut.max_tension] = StructuresAnalysis(lower_strut);

%% Tanks
[structural_loads.fuel_tank.max_compression, structural_loads.fuel_tank.max_tension] = StructuresAnalysis(FuelTank);
[structural_loads.oxygen_tank.max_compression, structural_loads.oxygen_tank.max_tension] = StructuresAnalysis(OxyTank);
[structural_loads.copv_tube.max_compression, structural_loads.copv_tube.max_tension] = StructuresAnalysis(COPVTube);

LBF2N = 4.44822;

% convert lbf values to newtons
component_names = fieldnames(structural_loads);
for i = 1:length(component_names)
    component_name = component_names{i};
    load_names = fieldnames(structural_loads.(component_name));

    for j = 1:length(load_names)
        load_name = load_names{j};
        structural_loads.(component_name).(load_name) = structural_loads.(component_name).(load_name) * LBF2N;
    end
end

save("structural_loads.mat", "structural_loads");
