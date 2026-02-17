classdef Material_Properties_Class
    properties
        name (1,1) string;
        youngs_modulus (1,1) double; % (psi)
        yieldCompressionStrength (1,1) double; % (psi)
        yieldTensionStrength (1,1) double; % (psi) 
        ultimateTensionStrength (1,1) double; % (psi)
        density (1,1) double % (lb/in^3)
    end
end