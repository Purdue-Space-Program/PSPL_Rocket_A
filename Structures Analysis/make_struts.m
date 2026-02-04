obtain_vehicle_parameters
parameters = load("vehicle_parameters.mat").vehicle_parameters;
%%
lower_strut = StrutClass;
lower_strut.name = "Lower Strut";
lower_strut.shape = 'Asym T';
lower_strut.width = 3/4; % [in]
lower_strut.wallThickness = 1/8; % [in]
lower_strut.length = parameters.lower_length;
lower_strut.distance = 22.6535; % [in]
lower_strut.radius = 2; % [in]
lower_strut.material = Aluminum6061T6MaterialProperties;
% from cad
lower_strut.crossArea = 0.5837; % [in^2]
lower_strut.radiusGyration = 2.6247; % [in]
%%
mid_strut = StrutClass;
mid_strut.name = 'Mid Strut';
mid_strut.shape = 'Circle';
mid_strut.oD = 3/4;
mid_strut.iD = 1/2;
mid_strut.length = 5; % Length of strut (in)
mid_strut.distance = 39.0935; % Location of strut from aft, top of strut (in)
mid_strut.radius = 2.25; % Distance from center axis (in) 
mid_strut.material = Aluminum6063T52MaterialProperties; % Material
%%
upper_strut = StrutClass;
upper_strut.name = 'Upper Strut';
upper_strut.shape = 'Circle';
upper_strut.oD = 3/4; % Width, of the side weak to bending (in)
upper_strut.iD = 1/8; % Wall thickness, of the side weak to bending (in)
upper_strut.length = 5.7; % Length of strut (in)
upper_strut.distance = 46.0935 + 5.7; % Location of strut from aft, top of strut (in)
upper_strut.radius = 2.25; % Distance from center axis (in) 
upper_strut.material = Aluminum6063T52MaterialProperties;



% classdef LowerStrutValues
%     properties
%         name = 'Lower Strut'
%         width = 3/4; % Width, of the side weak to bending (in)
%         wallThickness = 1/8; % Wall thickness, of the side weak to bending (in)
%         length = 14.45 % Length of strut (in)
%         distance = 22.6535 % Location of strut from aft, top of strut (in)
%         radius = 2; % Distance from center axis (in) 
%         crossArea = 0.5837; % Cross-sectional area pulled from NX (in^2)
%         radiusGyration = 2.6247; % Radius of gyration pulled from NX (in)
%         material = Aluminum6063T52MaterialProperties
%     end
% end