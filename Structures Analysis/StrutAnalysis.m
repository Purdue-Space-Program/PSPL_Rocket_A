%function StrutAnalysis(objectInQuestion)
%% ____________________
%% INITIALIZATION

%% SFD Importing
currentDirectory = pwd;
cd('C:../SFD')

sfdData = load('sfd_outputs_max_q.mat');
%sfdData = load('sfd_outputs_off_the_rail.mat');
axialLoads = sfdData.axial_array / 4.448; % Axial loads converted to pounds
shearLoads = sfdData.shear_array / 4.448; % Shear loads converted to pounds
momentLoads = sfdData.bending_array * 8.85; % Moment loads converted to inch-pounds
lengthLoads = sfdData.length_along_rocket_linspace * 39.37; % Length along the rocket converted inches

cd(currentDirectory)

%% Object Properities

objectInQuestion = MidStrutValues;

%% ____________________
%% Parameters
safetyFactor = [1.5, 1.5]; % Safety Factor [compression, tension] (randomly chosen fr fr)
K = 0.9; % Effective length factor Fixed-Fixed

length = objectInQuestion.length;
radiusGyration = objectInQuestion.radiusGyration;
crossArea = objectInQuestion.crossArea;
material = objectInQuestion.material;
distance = objectInQuestion.distance;
radius = objectInQuestion.radius;
name = objectInQuestion.name;

%% ____________________
%% Calculated Properities
effectiveLength = K * length;
slendernessRatio = effectiveLength / radiusGyration;
mass = length * crossArea * material.density;

buckleLimit = (material.yieldCompressionStrength - (material.yieldCompressionStrength * K * effectiveLength / (2 * pi * radiusGyration)) ^ 2 * (material.youngs ^ (-1))) * crossArea;
compressionLimit = material.yieldCompressionStrength * crossArea; % Compressive limit
tensionLimit = material.yieldTensionStrength * crossArea; % Tension limit
eulerLimit = (pi ^ 2 * material.youngs / ((K * effectiveLength / radiusGyration) ^ 2)) * crossArea;


%% SFD Properities
[~, location] = min(abs(lengthLoads - distance));
axialLoad = axialLoads(location);
torqueLoad = momentLoads(location);

%% Load Limit Properities
netLoad = [axialLoad + 2 * torqueLoad / radius, axialLoad - 2 * torqueLoad / radius] / 3; % Net force
failureMode = min([buckleLimit, eulerLimit, compressionLimit]);

safetyAllowance(1) = (failureMode - (safetyFactor(1) * netLoad(1))) / (safetyFactor(1) * netLoad(1));
safetyAllowance(2) = (tensionLimit - abs(safetyFactor(2) * netLoad(2))) / abs((safetyFactor(2) * netLoad(2)));
%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS
fprintf("Analysis of %s\n", name)
fprintf("Max compression case: %.2f lbf\n", netLoad(1));
if netLoad(2) >= 0
    fprintf("Max tension case: Not Under Tension\n\n");
else
    fprintf("Max tension case: %.2f lbf\n\n", netLoad(2));
end
fprintf("Slenderness Ratio: %.2f\n", slendernessRatio)
fprintf("Buckling Load Limit: %.2f lbf\n", buckleLimit)
%fprintf("Euler Buckling Load Limit: %.2f lbf\n", eulerLimit)
fprintf("Compression Load Limit: %.2f lbf\n", compressionLimit)
fprintf("Tension Load Limit: %.2f lbf\n", tensionLimit)
fprintf("Max compressive load safety factor for the %s: %.2f\n", name, safetyAllowance(1))
fprintf("Max tension load safety factor for the %s: %.2f\n", name, safetyAllowance(2))
fprintf("------------------------------------------------------\n")

%% ____________________
%% Functions
% function [radiusGyration, crossArea] = squareRadiusGyrationCalc(sideLength, thickness)
% shortSide = sideLength - 2*thickness;
% moment = (sideLength^4 - shortSide^4) / 12;
% crossArea = sideLength^2 - shortSide^2;
% radiusGyration = (moment / crossArea) ^ (1/2);
% end
% 
% function limit = bucklingEuler(effectiveLength, youngsMod, K, radiusGyration, crossArea)
% sigma = pi ^ 2 * youngsMod / ((K * effectiveLength / radiusGyration) ^ 2);
% limit = sigma * crossArea;
% end
% 
% function limit = bucklingJohnson(yieldStrength, effectiveLength, youngsMod, K, radiusGyration, crossArea)
% sigma = yieldStrength - (yieldStrength * K * effectiveLength / (2 * pi * radiusGyration)) ^ 2 * (youngsMod ^ (-1));
% limit = sigma * crossArea;
% end
% 
% function location = lengthLocation(target, lengthLoads)
% differences = abs(lengthLoads - target);
% [~, location] = min(differences);
% %closestValue = myArray(indexOfClosest);
% end