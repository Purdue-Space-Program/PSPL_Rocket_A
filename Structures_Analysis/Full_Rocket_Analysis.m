clear;
% clc;

structural_loads = struct();

make_components

%% Struts

[structural_loads.upper_strut.max_compression, structural_loads.upper_strut.max_tension] = StrutucturesAnalysis(upper_strut);
[structural_loads.mid_strut.max_compression, structural_loads.mid_strut.max_tension] = StrutucturesAnalysis(mid_strut);
[structural_loads.lower_strut.max_compression, structural_loads.lower_strut.max_tension] = StrutucturesAnalysis(lower_strut);

%% Tanks
[structural_loads.fuel_tank.max_compression, structural_loads.fuel_tank.max_tension] = StrutucturesAnalysis(FuelTank);
[structural_loads.oxygen_tank.max_compression, structural_loads.oxygen_tank.max_tension] = StrutucturesAnalysis(OxyTank);
[structural_loads.copv_tube.max_compression, structural_loads.copv_tube.max_tension] = StrutucturesAnalysis(COPVTube);

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
