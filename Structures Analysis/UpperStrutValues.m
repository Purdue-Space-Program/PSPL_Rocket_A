classdef UpperStrutValues
    properties
        name = 'Upper Strut'
        shape = 'Circle'
%        width = 3/4; % Width, of the side weak to bending (in)
%        wallThickness = 1/8; % Wall thickness, of the side weak to bending (in)
        oD = 3/4
        iD = 1/2
        length = 5.7 % Length of strut (in)
        distance = 46.0935 + 5.7 % Location of strut from aft, top of strut (in)
        radius = 2.25; % Distance from center axis (in) 
%        crossArea = (3/4)^2 - (3/4 - 2 * 1/8)^2; % (in)
%        radiusGyration = SquareRadiusGyrationCalc(3/4, 1/8); % (in)
        material = Aluminum6063T52MaterialProperities % Material
    end
end