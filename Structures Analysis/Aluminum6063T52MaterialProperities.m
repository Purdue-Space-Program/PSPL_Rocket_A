classdef Aluminum6063T52MaterialProperities
    properties
        youngs = 10007601 % (psi)
        yieldCompressionStrength = 31000; % (psi) check please
        yieldTensionStrength = 35000; % (psi) bad
        ultimateTensionStrength = 42000; % (psi) bad
        density = 0.0975 % (lb/in^3) bad
    end
end
% Sourced from pg 857 of mmpds