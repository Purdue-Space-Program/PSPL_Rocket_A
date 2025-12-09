classdef AluminumMaterialProperities
    properties
        youngs = 10007601 % (psi)
        yieldCompressionStrength = 35000; % (psi)
        yieldTensionStrength = 35000; % (psi) 
        ultimateTensionStrength = 42000; % (psi)
        density = 0.0975 % (lb/in^3)
    end
end
% Sourced from pg 857 of mmpds