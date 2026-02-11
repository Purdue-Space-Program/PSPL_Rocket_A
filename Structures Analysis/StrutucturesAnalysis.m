function [maxCompression, maxTension] = StrutucturesAnalysis(objectInQuestion)
%% ____________________
%% INITIALIZATION

%% Object Properties

%objectInQuestion = MidStrutValues;

%% ____________________
%% Parameters
safetyFactorTension = 1.5; % Honest guess at a safety factor
safetyFactorBending = 1.65; % Bending safety factor (ADM ASD)

K = 0.9; % Effective length factor Fixed-Fixed

length = objectInQuestion.length; % Length (in)
distance = [objectInQuestion.distance, objectInQuestion.distance - length]; % Distance from Aft to top (in)
radius = objectInQuestion.radius; % Distance from center axis (in)
name = objectInQuestion.name; % Name of the objectin in question
material = objectInQuestion.material; % Material properties from object


if objectInQuestion.shape == "Square"
    width = objectInQuestion.width; % Width of the object (in)
    wallThickness = objectInQuestion.wallThickness; % Wall thickness of the object (in)
    crossArea = width^2 - ((width - 2 * wallThickness)^2); % Cross sectional area)
    radiusGyration = SquareRadiusGyrationCalc(width, wallThickness); % Radius of gyration

elseif objectInQuestion.shape == "Circle"
    OD = objectInQuestion.OD;
    ID = objectInQuestion.ID;
    crossArea = ((OD/2)^2 - (ID/2)^2) * pi;
    radiusGyration = ((((OD)^4 - (ID)^4) / 64) / crossArea) ^ (1/2);
    width = OD; % Width of the object (in)
    wallThickness = (OD - ID)/2; % Wall thickness of the object (in)
else
    width = objectInQuestion.width; % Width of the object (in)
    wallThickness = objectInQuestion.wallThickness; % Wall thickness of the object (in)
    crossArea = objectInQuestion.crossArea; % Cross sectional area)
    radiusGyration = objectInQuestion.radiusGyration; % Radius of gyration
end
%% ____________________
%% Calculated Properties

effectiveLength = K * length;
slendernessRatio = effectiveLength / radiusGyration;
mass = length * crossArea * material.density;

transitionalRatio = sqrt(2 * pi^2 * material.youngs_modulus) / (K ^2 * material.yieldCompressionStrength);
wtRatio = width / (wallThickness); % Width-Thickness ratio

buckleLimit = (material.yieldCompressionStrength - (material.yieldCompressionStrength * K * effectiveLength / (2 * pi * radiusGyration)) ^ 2 * (material.youngs_modulus^ (-1))) * crossArea;
compressionLimit = material.yieldCompressionStrength * crossArea; % Compressive limit
tensionLimit = material.yieldTensionStrength * crossArea; % Tension limit
eulerLimit = (pi ^ 2 * material.youngs_modulus / ((K * effectiveLength / radiusGyration) ^ 2)) * crossArea;

%% Buckling Flow

if wtRatio <= ((0.11 * material.youngs_modulus) / material.yieldCompressionStrength)
    area = crossArea;
    Fe = (pi^2 * material.youngs_modulus) / (slendernessRatio ^ 2);
    consideringLocalBuckling = 1;
else
    area = crossArea * ((0.038 * material.youngs_modulus) / (material.yieldCompressionStrength * (wtRatio)) + 2/3);
    Fe = (pi ^ 2 * material.youngs_modulus) / (slendernessRatio ^ 2);
    consideringLocalBuckling = 0;
end

if slendernessRatio <= (4.71 * (material.youngs_modulus / material.yieldCompressionStrength)^ (1/2))
    Fcr = ((0.658) ^ (material.yieldCompressionStrength / Fe)) * material.yieldCompressionStrength;
else
    Fcr = 0.877 * Fe;
end

pAllow = (Fcr * area) / safetyFactorBending; % Allowable axial load [lb]

%% Load Limit Properties

[maxCompression, maxTension] = NetAxialLoad(distance, radius); % Max compressive and tensile cases, divided by the three struts
maxCompression = maxCompression / 3;
maxTension = abs(maxTension) / 3;

% === Tension Check ===
tensileStress = maxTension / crossArea;
MoSCompression = (pAllow / maxCompression) - 1;
MoSTension = ((area * material.yieldTensionStrength) / (tensileStress * safetyFactorTension)) - 1;


%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS
fprintf("Analysis of %s\n", name)
fprintf("\tMaterial: %s\n", objectInQuestion.material.name)
fprintf("\tShape: %s\n", objectInQuestion.shape)
fprintf("\tWidth: %.2f\n", width)
fprintf("\tWall Thickness: %.2f\n", wallThickness)
fprintf("\tMass of %s: %.2f lbs\n", name, mass);

fprintf("\n\tAvailable Axial Strength (ASD): %.2f lbs\n", pAllow);
fprintf("\tMax compression case: %.2f lbf\n", maxCompression);
fprintf("\tMax tension case: %.2f lbf\n", maxTension);

if consideringLocalBuckling == 1
    fprintf("\tNon-slender - local buckling not critical.\n")
else
    fprintf("\tSlender = local buckling must be considered.\n");
end
% fprintf("Slenderness Ratio: %.2f\n", slendernessRatio)
% fprintf("Buckling Load Limit: %.2f lbf\n", buckleLimit)
% %fprintf("Euler Buckling Load Limit: %.2f lbf\n", eulerLimit)
% fprintf("Compression Load Limit: %.2f lbf\n", compressionLimit)
% fprintf("Tension Load Limit: %.2f lbf\n", tensionLimit)
fprintf("\tCompression MoS: %.2f\n", MoSCompression);
fprintf("\tTension Yield MoS: %.2f\n", MoSTension);
fprintf("------------------------------------------------------\n\n\n")
% fprintf("Max compressive load safety factor for the %s: %.2f\n", name, safetyAllowance(1))
% fprintf("Max tension load safety factor for the %s: %.2f\n", name, safetyAllowance(2))
% fprintf("------------------------------------------------------\n")

%% ____________________
%% Functions
% function [radiusGyration, crossArea] = squareRadiusGyrationCalc(sideLength, thickness)
% shortSide = sideLength - 2*thickness;
% moment = (sideLength^4 - shortSide^4) / 12;
% crossArea = sideLength^2 - shortSide^2;
% radiusGyration = (moment / crossArea) ^ (1/2);
% end
% 
% function limit = bucklingEuler(effectiveLength, youngs_modulus, K, radiusGyration, crossArea)
% sigma = pi ^ 2 * youngs_modulus / ((K * effectiveLength / radiusGyration) ^ 2);
% limit = sigma * crossArea;
% end
% 
% function limit = bucklingJohnson(yieldStrength, effectiveLength, youngs_modulus, K, radiusGyration, crossArea)
% sigma = yieldStrength - (yieldStrength * K * effectiveLength / (2 * pi * radiusGyration)) ^ 2 * (youngs_modulus ^ (-1));
% limit = sigma * crossArea;
% end
% 
% function location = lengthLocation(target, lengthLoads)
% differences = abs(lengthLoads - target);
% [~, location] = min(differences);
% %closestValue = myArray(indexOfClosest);
% end