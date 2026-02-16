clear;
clc;

structural_loads = struct();

%% Struts
make_struts

%[structural_loads.upper_strut.max_compression, structural_loads.upper_strut.max_tension] = StrutucturesAnalysis(upper_strut);
%[structural_loads.mid_strut.max_compression, structural_loads.mid_strut.max_tension] = StrutucturesAnalysis(mid_strut);
[structural_loads.lower_strut.max_compression, structural_loads.lower_strut.max_tension] = StrutucturesAnalysis(lower_strut);

%% Tanks
%make_tanks

%[structural_loads.fuel_tank.max_compression, structural_loads.fuel_tank.max_tension] = StrutucturesAnalysis(FuelTank);
%[structural_loads.oxygen_tank.max_compression, structural_loads.oxygen_tank.max_tension] = StrutucturesAnalysis(OxyTank);
%[structural_loads.copv_tank.max_compression, structural_loads.copv_tank.max_tension] = StrutucturesAnalysis(COPVWall);

%save("structural_loads.mat", "structural_loads");
