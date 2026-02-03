classdef MidStrutValues
    properties
        name = 'Mid Strut'
        shape = 'Circle'
%        width = 3/4; % Width, of the side weak to bending (in)
%        wallThickness = 1/8; % Wall thickness, of the side weak to bending (in)
        oD = 3/4
        iD = 1/2
        length = 5 % Length of strut (in)
        distance = 39.0935 % Location of strut from aft, top of strut (in)
        radius = 2.25; % Distance from center axis (in) 
%        crossArea = ((((3/4)/2)^2 - ((1/2)/2)^2) * pi); % (in)
%        radiusGyration = ((((3/4)^4 - (1/2)^4) / 64) / ((((3/4)/2)^2 - ((1/2)/2)^2) * pi)) ^ (1/2); % (in)
        material = Aluminum6063T52MaterialProperities % Material
    end
end