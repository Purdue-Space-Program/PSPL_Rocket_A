classdef OxidizerTankValues
    properties
        name = 'Oxidizer Tank'
        length = 8 % Length of strut (in)
        distance = 22.6535 % Location of strut from aft, top of strut (in)
        radius = 3; % Distance from center axis (in) 
        iD = 5.75; % (in)
        oD = 6; % (in)
        material = Aluminum6061T6MaterialProperties
    end
end