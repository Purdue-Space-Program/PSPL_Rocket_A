make_material_properties
obtain_vehicle_parameters
parameters = load("vehicle_parameters.mat").vehicle_parameters;
wet_mass_distribution = load("wet_mass_distribution.mat").wet_mass_distribution;


IN2M = 0.0254;  % [m/in] Conversion factor from in to m
M2IN = 1 / IN2M;  % [in/m] Conversion factor from m to in

M2FT = 3.28084;  % [ft/m] Conversion factor from m to ft
FT2M = 1 / M2FT;  % [m/ft] Conversion factor from ft to m


FuelTank = TankClass;
FuelTank.name = 'Fuel Tank';
FuelTank.shape = "Circle";
FuelTank.OD = 6; % [in]
FuelTank.ID = 5.75; % [in]
FuelTank.length = 8; % from CAD (in)
FuelTank.distance = 43.6535; % Location of strut from aft, top of strut (in)
FuelTank.radius = 3; % Distance from center axis (in) 
FuelTank.material = Aluminum_6061_T6_Material_Properties;


OxyTank = TankClass;
OxyTank.name = 'Oxidizer Tank';
OxyTank.shape = "Circle";
OxyTank.OD = 6; % [in]
OxyTank.ID = 5.75; % [in]
OxyTank.length = 8; % from CAD (in)
OxyTank.distance = 22.6535; % Location of strut from aft, top of strut (in)
OxyTank.radius = 3; % Distance from center axis (in) 
OxyTank.material = Aluminum_6061_T6_Material_Properties;


COPVWall.name = 'COPV Wall';
COPVWall.shape = "Circle";
COPVWall.OD = 6; % [in]
COPVWall.ID = 5.75; % [in]
COPVWall.length = 5; % from CAD (in)
COPVWall.distance = 40; % Location of strut from aft, top of strut (in)
COPVWall.radius = 3; % Distance from center axis (in) 
COPVWall.material = Aluminum_6061_T6_Material_Properties;

