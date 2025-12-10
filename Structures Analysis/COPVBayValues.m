classdef COPVBayValues
    properties
        name = 'COPV Tank'
        length = 5 % Length of strut (in)
        distance = 40 % Location of strut from aft, top of strut (in)
        radius = 0; % Distance from center axis (in) 
        effectiveLength % Effective length of the strut
        radiusGyration = SquareRadiusGyrationCalc(0.5, 0.065);
        crossArea = (0.5)^2 - (0.5 - 0.065*2)^2;
        slendernessRatio % Slenderness Ratio
        buckleLimit % Buckling Limit (lbf)
        location % Location in the SFD calculations
        axialLoad % Axials load (lbf)
        torqueLoad % Moment load (lb-in)
        netLoad % Net load on strut (lbf)
        material = AluminumMaterialProperities % Material of the
        mass % Mass of the strut
        safetyAllowance % Amount of safety factor
    end
end