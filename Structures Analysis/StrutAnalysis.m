function StrutAnalysis(objectInQuestion)
%% ____________________
%% INITIALIZATION

%% SFD Importing
currentDirectory = pwd;
cd('C:../SFD')

sfdData = load('sfd_outputs_max_q.mat');
%sfdData = load('sfd_outputs_off_the_rail.mat');
axialLoads = sfdData.axial_array / 4.448; % Axial loads converted to pounds
shearLoads = sfdData.shear_array / 4.448; % Shear loads converted to pounds
momentLoads = sfdData.bending_array * 8.85; % Moment loads converted to feet-pounds
lengthLoads = sfdData.length_along_rocket_linspace * 39.37; % Length along the rocket converted inches

cd(currentDirectory)

%% Object Properities

%% ____________________
%% Parameters
safetyFactor = 1.4; % Safety Factor (randomly chosen fr fr)
K = 0.9; % Effective length factor Fixed-Fixed

%% ____________________
%% Calculated Properities
objectInQuestion.effectiveLength = K * objectInQuestion.length;
objectInQuestion.slendernessRatio = objectInQuestion.effectiveLength / objectInQuestion.radiusGyration;
objectInQuestion.mass = objectInQuestion.length * objectInQuestion.crossArea * objectInQuestion.material.density;

objectInQuestion.buckleLimit = (objectInQuestion.material.yieldStrength - (objectInQuestion.material.yieldStrength * K * objectInQuestion.effectiveLength / (2 * pi * objectInQuestion.radiusGyration)) ^ 2 * (objectInQuestion.material.youngs ^ (-1))) * objectInQuestion.crossArea;

%% SFD Properities
[~, objectInQuestion.location] = min(abs(lengthLoads - objectInQuestion.distance));
objectInQuestion.axialLoad = axialLoads(objectInQuestion.location);
objectInQuestion.torqueLoad = momentLoads(objectInQuestion.location);

%% Load Limit Properities
objectInQuestion.netLoad = [objectInQuestion.axialLoad + 2 * objectInQuestion.torqueLoad / objectInQuestion.radius, objectInQuestion.axialLoad - 2 * objectInQuestion.torqueLoad / objectInQuestion.radius]; % Net force
objectInQuestion.safetyAllowance = objectInQuestion.buckleLimit ./ (safetyFactor * objectInQuestion.netLoad(1));

%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS
fprintf("Analysis of %s\n", objectInQuestion.name)
fprintf("Max compression case: %.2f lbf\n", objectInQuestion.netLoad(1));
if objectInQuestion.netLoad(2) >= 0
    fprintf("Max tension case: Not Under Tension\n\n");
else
    fprintf("Max tension case: %.2f lbf\n\n", objectInQuestion.netLoad(2));
end
fprintf("Buckling Load Limit: %.2f lbf\n", objectInQuestion.buckleLimit)
fprintf("Max load safety factor for the %s: %.2f\n", objectInQuestion.name, objectInQuestion.safetyAllowance - 1)
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