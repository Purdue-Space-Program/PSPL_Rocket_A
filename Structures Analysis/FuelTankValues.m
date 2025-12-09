classdef FuelTankValues
    properties
        name = 'Fuel tank'
        length = 8 % Length of strut (in)
        distance = 43.6535 % Location of strut from aft, top of strut (in)
        radius = 0; % Distance from center axis (in) 
        effectiveLength % Effective length of the strut
        radiusGyration = ((3^2 + (3-1/8)^2) / 2)^(1/2);
        crossArea = ((3)^2 - (3 - 1/8)^2) * pi;
        slendernessRatio % Slenderness Ratio
        buckleLimit % Buckling Limit (lbf)
        location % Location in the SFD calculations
        axialLoad % Axials load (lbf)
        torqueLoad % Moment load (lb-in)
        netLoad % Net load on strut (lbf)
        material = SteelMaterialProperities % Material of the
        mass % Mass of the strut
        safetyAllowance % Amount of safety factor
    end
end