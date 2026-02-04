classdef StrutClass
    properties
        name (1,1) string;
        width (1,1) double; % Width, of the side weak to bending (in)
        wallThickness (1,1) double; % Wall thickness, of the side weak to bending (in)
        length (1,1) double % Length of strut (in)
        distance (1,1) double % Location of strut from aft, top of strut (in)
        radius (1,1) double; % Distance from center axis (in) 
        crossArea (1,1) double; % (in)
        radiusGyration (1,1) double; % (in)
        material = Aluminum6061T6MaterialProperities % Material
    end
end