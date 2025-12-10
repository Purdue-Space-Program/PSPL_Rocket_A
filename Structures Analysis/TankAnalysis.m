%function StrutAnalysis(objectInQuestion)
%% ____________________
%% INITIALIZATION
clear
clc

%% SFD Importing
currentDirectory = pwd;
cd('C:../SFD')

sfdData = load('sfd_outputs_max_q.mat');
%sfdData = load('sfd_outputs_off_the_rail.mat');
axialLoadssfd = sfdData.axial_array / 4.448; % Axial loads converted to pounds
shearLoadssfd = sfdData.shear_array / 4.448; % Shear loads converted to pounds
momentLoadssfd = sfdData.bending_array * 8.85; % Moment loads converted to inch-pounds
lengthLoadssfd = sfdData.length_along_rocket_linspace * 39.37; % Length along the rocket converted inches

rfdData = load('sfd_outputs_max_q.mat');
axialLoadsrfd = rfdData.axial_array / 4.448: % Axial loads converted to pounds
shearLoadsrfd = rfdData.shear_array / 4.448; % Shear loads converted to pounds
momentLoadsrfd = rfdData.bending_array * 8.85; % Moment loads conevrted to inch-pounds
lengthLoadsrfd = rfdData.length_along_rocket_linspace * 39.37; % Length along the rocket converted inches

cd(currentDirectory)

%% Object Properities

objectInQuestion = FuelTankValues;

%% ____________________
%% Parameters
safetyFactor = [1.5, 1.5]; % Safety Factor [compression, tension] (randomly chosen fr fr)
safetyFactorComp = 1.65; % Compression safety factor (ADM ASD)
safetyFactorBending = 1.67; % Bending safety factor (ADM ASD)

K = 0.9; % Effective length factor Fixed-Fixed

iD = objectInQuestion.iD; % Inner diameter (in)
oD = objectInQuestion.oD; % Outer diameter (in)
length = objectInQuestion.length; % Length (in)
distance = objectInQuestion.distance; % Distance from Aft to top (in)
radius = objectInQuestion.radius; % Distance from center axis (in)
name = objectInQuestion.name; % Name of the objectin in question

material = objectInQuestion.material; % Material properties from object


%% ____________________
%% Calculated Properities
crossArea = ((oD/2)^2 - (iD/2)^2) * pi;
radiusGyration = ((((oD)^4 - (iD)^4) / 64) / crossArea) ^ (1/2);

effectiveLength = K * length;
slendernessRatio = effectiveLength / radiusGyration;
mass = length * crossArea * material.density;

transitionalRatio = sqrt(2 * pi^2 * material.youngs) / (K ^2 * material.yieldCompressionStrength);
wtRatio = oD / (oD - iD); % Width-Thickness ratio
plasticMod = (pi/4) * ((oD/2) ^ 3 - (iD/2)^3); % Plastic modulus
sectionMod = (pi * (oD^4 - iD^4)) / (32 * oD); % Section modulus

buckleLimit = (material.yieldCompressionStrength - (material.yieldCompressionStrength * K * effectiveLength / (2 * pi * radiusGyration)) ^ 2 * (material.youngs ^ (-1))) * crossArea;
compressionLimit = material.yieldCompressionStrength * crossArea; % Compressive limit
tensionLimit = material.yieldTensionStrength * crossArea; % Tension limit
eulerLimit = (pi ^ 2 * material.youngs / ((K * effectiveLength / radiusGyration) ^ 2)) * crossArea;

%% Buckling Constants?
if wtRatio <= ((0.11 * material.youngs) / material.yieldCompressionStrength)
    area = crossArea;
    Fe = (pi^2 * material.youngs) / (slendernessRatio ^ 2);
    consideringLocalBuckling = 1;
else
    area = crossArea * ((0.038 * material.youngs) / (material.yieldCompressionStrength * (wt_ratio)) + 2/3);
    Fe = (pi ^ 2 * E) / (slenderness_ratio ^ 2);
    consideringLocalBuckling = 0;
end

if slendernessRatio <= (4.71 * (material.youngs / material.yieldCompressionStrength)^ (1/2))
    Fcr = ((0.658) ^ (material.yieldCompressionStrength / Fe)) * material.yieldCompressionStrength;
else
    Fcr = 0.877 * Fe;
end

pAllow = (Fcr * area) / safetyFactorComp; % Allowable axial load (lb)

% === Combined Axial and Bending Check ===
area_total = (pi/4) * (oD^2 - iD^2);

%% SFD Properities
[~, locationTakeOff] = min(abs(lengthLoadssfd - distance));
[~, locationRecovery] = min(abs(lengthLoadss - distance));

%% Load Limit Properities
netLoadTakeoff = [axialLoadssfd(locationTakeOff) + 2 * momentLoadssfd(locationTakeOff) / radius, axialLoadssfd(locationTakeOff) - 2 * momentLoadssfd(locationTakeOff) / radius]; % Net force
netLoadRecovery = [];

netLoad(1) = max([]);
netLoad(2) = max([]);

% === Tension Check ===
tensileStress = abs(netLoad(2) / area_total);
FoSTension = material.yieldCompressionStrength / tensileStress - 1;

%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS
fprintf("Analysis of %s\n", name)
fprintf("Max compression case: %.2f lbf\n", netLoad(1));
fprintf("Max tension case: %.2f lbf\n\n", netLoad(2));

fprintf("Strut OD: %.2f in, ID: %.2f in\n", oD, iD);
if consideringLocalBuckling == 1
    fprintf("Nonslender - local buckling not critcal.\n")
else
    fprintf("Slender = local buckling must be considered.\n");
end
% fprintf("Slenderness Ratio: %.2f\n", slendernessRatio)
% fprintf("Buckling Load Limit: %.2f lbf\n", buckleLimit)
% %fprintf("Euler Buckling Load Limit: %.2f lbf\n", eulerLimit)
% fprintf("Compression Load Limit: %.2f lbf\n", compressionLimit)
% fprintf("Tension Load Limit: %.2f lbf\n", tensionLimit)
fprintf("Available Axial Strength (ASD): %.2f lbs\n", pAllow);
fprintf("Mass of %s: %.2f lbs\n", name, mass);
fprintf("Tensile FoS: %.2f\n", FoSTension);
fprintf("------------------------------------------------------\n")
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