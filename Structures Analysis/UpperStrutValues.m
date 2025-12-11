classdef UpperStrutValues
    properties
        name = 'Upper Strut'
        width = 3/4; % Width, of the side weak to bending (in)
        wallThickness = 1/8; % Wall thickness, of the side weak to bending (in)
        length = 5 % Length of strut (in)
        distance = 48.7111 + 4 % Location of strut from aft, top of strut (in)
        radius = 2.25; % Distance from center axis (in) 
        crossArea = (3/4)^2 - (3/4 - 2 * 1/8)^2; % (in)
        radiusGyration = SquareRadiusGyrationCalc(3/4, 1/8); % (in)
        material = Aluminum6063T52MaterialProperities % Material
    end
end
% This strut isn't finished as of 12/7 at 11:53 am, 