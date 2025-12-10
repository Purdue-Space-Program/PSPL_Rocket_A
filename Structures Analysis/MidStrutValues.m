classdef MidStrutValues
    properties
        name = 'Mid Strut'
        length = 5 % Length of strut (in)
        distance = 35.6535 % Location of strut from aft, top of strut (in)
        radius = 2.25; % Distance from center axis (in) 
        crossArea = .12; % (in)
        radiusGyration = 2.8404; % (in)
        material = Aluminum6063T52MaterialProperities % Material of the
    end
end