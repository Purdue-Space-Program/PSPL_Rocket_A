function [maxCompression, maxTension] = NetAxialLoad(location, radius, graphStatus)
%
%
% THIS FUNCTION OUTPUTS THE FORCES RAW, YOU STILL NEED TO DIVIDE BY 3 FOR
% STRUTS. This function works for any point on the rocket.This function
% returns the maximum force that could appear at that location according to
% the latest running of the SFD and RFD. 


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
axialLoadsrfd = rfdData.axial_array_recovery / 4.448; % Axial loads converted to pounds
shearLoadsrfd = rfdData.shear_array_recovery / 4.448; % Shear loads converted to pounds
momentLoadsrfd = rfdData.bending_array_recovery * 8.85; % Moment loads conevrted to inch-pounds
lengthLoadsrfd = rfdData.length_along_rocket_linspace * 39.37; % Length along the rocket converted inches

cd(currentDirectory)

%% Net Load Calcs

topLocation = max(location);
bottomLocation = min(location);

[~, locationTakeOff(1)] = min(abs(lengthLoadssfd - topLocation));
[~, locationTakeOff(2)] = min(abs(lengthLoadssfd - bottomLocation));
[~, locationRecovery(1)] = min(abs(lengthLoadsrfd - topLocation));
[~, locationRecovery(2)] = min(abs(lengthLoadsrfd - bottomLocation));

if radius == 0
    netCompressionTakeoff = axialLoadssfd; 
    netTensionTakeoff = axialLoadssfd; 
    netCompressionRecovery = axialLoadsrfd; 
    netTensionRecovery = axialLoadsrfd;
else
netCompressionTakeoff = axialLoadssfd + 2 .* momentLoadssfd ./ radius; 
netTensionTakeoff = axialLoadssfd - 2 .* momentLoadssfd / radius; 
netCompressionRecovery = axialLoadsrfd + 2 .* momentLoadsrfd ./ radius; 
netTensionRecovery = axialLoadsrfd - 2 .* momentLoadsrfd / radius;
end

maxCompression = max([netCompressionTakeoff(locationTakeOff(1)), netCompressionRecovery(locationRecovery(1)), netCompressionTakeoff(locationTakeOff(2)), netCompressionRecovery(locationRecovery(2))]);
maxTension = min([netTensionTakeoff(locationTakeOff(1)), netTensionRecovery(locationRecovery(1)), netTensionTakeoff(locationTakeOff(2)), netTensionRecovery(locationRecovery(2))]);

%% GRAPHS

timeToGraph = exist('graphStatus');

if timeToGraph
    figure(1)
    subplot(1,3,1)
    plot(lengthLoadssfd, axialLoadssfd);
    hold on;
    plot([lengthLoadssfd(locationTakeOff(1)), lengthLoadssfd(locationTakeOff(2))], [axialLoadssfd(locationTakeOff(1)), axialLoadssfd(locationTakeOff(2))], '*r')
    hold off;
    grid on;
    title('Axial Loads at Takeoff');
    xlabel('Length along Rocket (inches)');
    ylabel('Axial Load (pounds)');

    subplot(1,3,2)
    plot(lengthLoadssfd, shearLoadssfd);
    hold on;
    plot([lengthLoadssfd(locationTakeOff(1)), lengthLoadssfd(locationTakeOff(2))], [shearLoadssfd(locationTakeOff(1)), shearLoadssfd(locationTakeOff(2))], '*r')
    hold off;
    grid on;
    title('Shear Loads at Takeoff');
    xlabel('Length along Rocket (inches)');
    ylabel('Shear Load (pounds)');

    subplot(1,3,3)
    plot(lengthLoadssfd, momentLoadssfd);
    hold on;
    plot([lengthLoadssfd(locationTakeOff(1)), lengthLoadssfd(locationTakeOff(2))], [momentLoadssfd(locationTakeOff(1)), momentLoadssfd(locationTakeOff(2))], '*r')
    hold off;
    grid on;
    title('Moment Loads at Takeoff');
    xlabel('Length along Rocket (inches)');
    ylabel('Moment Load (pounds-in)');


    figure(2)
    subplot(1,3,1)
    plot(lengthLoadsrfd, axialLoadsrfd);
    hold on;
    plot([lengthLoadsrfd(locationRecovery(1)), lengthLoadsrfd(locationRecovery(2))], [axialLoadsrfd(locationTakeOff(1)), axialLoadsrfd(locationTakeOff(2))], '*r')
    hold off;
    grid on;
    title('Axial Loads during Recovery');
    xlabel('Length along Rocket (inches)');
    ylabel('Axial Load (pounds)');

    subplot(1,3,2)
    plot(lengthLoadsrfd, shearLoadsrfd);
    hold on;
    plot([lengthLoadsrfd(locationRecovery(1)), lengthLoadsrfd(locationRecovery(2))], [shearLoadsrfd(locationTakeOff(1)), shearLoadsrfd(locationTakeOff(2))], '*r')
    hold off;
    grid on;
    title('Shear Loads at Recovery');
    xlabel('Length along Rocket (inches)');
    ylabel('Shear Load (pounds)');

    subplot(1,3,3)
    plot(lengthLoadsrfd, momentLoadsrfd);
    hold on;
    plot([lengthLoadsrfd(locationRecovery(1)), lengthLoadsrfd(locationRecovery(2))], [momentLoadsrfd(locationTakeOff(1)), momentLoadsrfd(locationTakeOff(2))], '*r')
    hold off;
    grid on;
    title('Moment Loads at Recovery');
    xlabel('Length along Rocket (inches)');
    ylabel('Moment Load (pounds-in)');

end