clear;
clc;

make_struts

%% Struts
[upper_strut_max_compression, upper_strut_max_tension] = StrutAnalysis(upper_strut);
[mid_strut_max_compression, mid_strut_max_tension] = StrutAnalysis(mid_strut);
[lower_strut_max_compression, lower_strut_max_tension] = StrutAnalysis(lower_strut);

%% Tanks
% TankAnalysis(FuelTankValues);
% TankAnalysis(OxidizerTankValues);
% add for copv tube?