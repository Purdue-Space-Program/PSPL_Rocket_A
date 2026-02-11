clear;
clc;

%% Struts
make_struts

[upper_strut_max_compression, upper_strut_max_tension] = StrutucturesAnalysis(upper_strut);
[mid_strut_max_compression, mid_strut_max_tension] = StrutucturesAnalysis(mid_strut);
[lower_strut_max_compression, lower_strut_max_tension] = StrutucturesAnalysis(lower_strut);

%% Tanks
make_tanks

[fuel_tank_max_compression, fuel_tank_max_tension] = StrutucturesAnalysis(FuelTank);
[oxy_tank_max_compression, oxy_tank_max_tension] = StrutucturesAnalysis(OxyTank);
[copv_tank_max_compression, copv_tank_max_tension] = StrutucturesAnalysis(COPVTank);