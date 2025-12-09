classdef SteelMaterialProperities
    properties
        youngs = 29000000 % (psi)
        yieldCompressionStrength = 33000; % (psi) % Currently a made up number
        yieldTensionStrength = 33000; % (psi) 
        ultimateTensionStrength = 45000; % (psi)
        density = 0.2818 % (lb/in^3)
    end
end
% Properities for grade C Cold-Formed Welded and Seamless Carbon Steel
% Structural Tubing in Rounds and Shapes
% https://compass.astm.org/content-access?contentCode=ASTM%7CA0500-03A%7Cen-US