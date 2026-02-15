make_material_properties
obtain_vehicle_parameters
parameters = load("vehicle_parameters.mat").vehicle_parameters;
wet_mass_distribution = load("wet_mass_distribution.mat").wet_mass_distribution;


IN2M = 0.0254;  % [m/in] Conversion factor from in to m
M2IN = 1 / IN2M;  % [in/m] Conversion factor from m to in

M2FT = 3.28084;  % [ft/m] Conversion factor from m to ft
FT2M = 1 / M2FT;  % [m/ft] Conversion factor from ft to m


upper_strut = StrutClass;
upper_strut.name = 'Upper Strut';
upper_strut.shape = "Square";
% upper_strut.OD = 3/4; % [in]
% upper_strut.ID = 1/2; % [in]
upper_strut.width = 3/4; % Width, of the side weak to bending [in]
upper_strut.wallThickness = 0.125; % Wall thickness, of the side weak to bending [in]
upper_strut.length = 11; % from CAD because parameters.upper_length isn't accurate since it goes to top of COPV tube [in]
upper_strut.distance = wet_mass_distribution.upper.top_distance_from_aft * M2IN; % Location of strut from aft to top of strut [in]
upper_strut.radius = 2.25; % Distance from center axis [in] 
upper_strut.material = Aluminum_6063_T52_Material_Properties;


mid_strut = StrutClass;
mid_strut.name = 'Mid Strut';
mid_strut.shape = "Square";
mid_strut.width = 3/4; % Width, of the side weak to bending [in]
mid_strut.wallThickness = 0.125; % Wall thickness, of the side weak to bending [in]
% mid_strut.OD = 1;
% mid_strut.ID = 1/2;
mid_strut.length = parameters.mid_length * M2IN; % Length of strut [in]
mid_strut.distance = wet_mass_distribution.mid.top_distance_from_aft * M2IN; % Location of strut from aft to top of strut [in]
mid_strut.radius = 2.25; % Distance from center axis [in] 
mid_strut.material = Aluminum_6063_T52_Material_Properties;


lower_strut = StrutClass;
lower_strut.name = "Lower Strut";
lower_strut.shape = "Asym T";
lower_strut.length = parameters.lower_length * M2IN;
lower_strut.distance = wet_mass_distribution.lower.top_distance_from_aft * M2IN; % Location of strut from aft to top of strut [in]
lower_strut.material = Aluminum_6061_T6_Material_Properties;

% from cad
lower_strut.crossArea = 0.49; % [in^2]
lower_strut.radiusGyration = 2.28239244; % [in]
lower_strut.radius = 1.962; % [in]
lower_strut.wallThickness = 0.2; % [in]
lower_strut.width = 1.65; % [in]
