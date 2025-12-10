classdef FuelTankValues
    properties
        name = 'Fuel tank'
        length = 8 % Length of strut (in)
        distance = 43.6535 % Location of strut from aft, top of strut (in)
        radius = 3; % Distance from center axis (in) 
        iD = 5.75; % (in)
        oD = 6; % (in)
        effectiveLength % Effective length of the strut
        crossArea % Cross sectional area (in^2)
        radiusGyration % Radius of gyration (idk)
        slendernessRatio % Slenderness Ratio
        buckleLimit % Buckling Limit (lbf)
        location % Location in the SFD calculations
        axialLoad % Axials load (lbf)
        torqueLoad % Moment load (lb-in)
        netLoad % Net load on strut (lbf)
        material = Aluminum6061T6MaterialProperities % Material of the
        mass % Mass of the strut
        safetyAllowance % Amount of safety factor
    end
end