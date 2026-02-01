classdef LowerStrutValues
    properties

        parameters = load("vehicle_parameters.mat");
        name = 'Lower Strut'
        width = 3/4; % Width, of the side weak to bending (in)
        wallThickness = 1/8; % Wall thickness, of the side weak to bending (in)
        length = 14.45; % Length of strut (in) 14.45 
        distance = 22.6535 % Location of strut from aft, top of strut (in)
        radius = 2; % Distance from center axis (in) 
        crossArea = 0.5837; % Cross-sectional area pulled from NX (in^2)
        radiusGyration = 2.6247; % Radius of gyration pulled from NX (in)
        material = Aluminum6063T52MaterialProperities % Material
    end
end