parameters = load("vehicle_parameters.mat");

lower_strut = StrutClass;

lower_strut.name = "Lower Strut";
lower_strut.width = 3/4; % [in]
lower_strut.wallThickness = 1/8; % [in]
lower_strut.length = parameters.lower_length;
lower_strut.distance = 22.6535; % [in]

% from cad
lower_strut.crossArea = 0.5837; % [in^2]
lower_strut.radiusGyration = 2.6247; % [in]

lower_strut.material = Aluminum6061T6MaterialProperities;

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
%         material = Aluminum6063T52MaterialProperities % Material
%     end
% end