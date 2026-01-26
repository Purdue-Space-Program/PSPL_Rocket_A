function [maxCompression, maxTension] = NetAxialLoad(location, radius)
%
%
%
%
%
%
% THIS FUNCTION OUTPUTS THE FORCES RAW, YOU STILL NEED TO DIVIDE BY 3 FOR
% STRUTS. This function works for any point on the rocket.This function
% returns the maximum force that could appear at that location according to
% the latest running of the SFD and RFD. 
%


%% Importing SFD and RFD
currentDirectory = pwd;
cd('C:../SFD')

sfdData = load('sfd_outputs_max_q.mat');
%sfdData = load('sfd_outputs_off_the_rail.mat');
axialLoadssfd = sfdData.axial_array / 4.448; % Axial loads converted to pounds
shearLoadssfd = sfdData.shear_array / 4.448; % Shear loads converted to pounds
momentLoadssfd = sfdData.bending_array * 8.85; % Moment loads converted to inch-pounds
lengthLoadssfd = sfdData.length_along_rocket_linspace * 39.37; % Length along the rocket converted inches

rfdData = load('rfd_outputs_recovery.mat');
axialLoadsrfd = rfdData.axial_array_45_deg / 4.448; % Axial loads converted to pounds
shearLoadsrfd = rfdData.shear_array_45_deg / 4.448; % Shear loads converted to pounds
momentLoadsrfd = rfdData.bending_array_45_deg * 8.85; % Moment loads conevrted to inch-pounds
lengthLoadsrfd = rfdData.length_along_rocket_linspace * 39.37; % Length along the rocket converted inches

cd(currentDirectory)

%% Net Load Calcs

[~, locationTakeOff] = min(abs(lengthLoadssfd - location));
[~, locationRecovery] = min(abs(lengthLoadsrfd - location));

netCompressionTakeoff = axialLoadssfd + 2 .* momentLoadssfd ./ radius; 
netTensionTakeoff = axialLoadssfd - 2 .* momentLoadssfd / radius; 
netCompressionRecovery = axialLoadsrfd + 2 .* momentLoadsrfd ./ radius; 
netTensionRecovery = axialLoadsrfd - 2 .* momentLoadsrfd / radius;

maxCompression = max([netCompressionTakeoff(locationTakeOff), netCompressionRecovery(locationRecovery)]);
maxTension = min([netTensionTakeoff(locationTakeOff), netTensionRecovery(locationRecovery)]);