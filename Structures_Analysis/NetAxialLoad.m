function [maxCompression, maxTension] = NetAxialLoad(location, radius, graphStatus)
%
%
% THIS FUNCTION OUTPUTS THE FORCES RAW, YOU STILL NEED TO DIVIDE BY 3 FOR
% STRUTS. This function works for any point on the rocket.This function
% returns the maximum force that could appear at that location according to
% the latest running of the SFD and RFD. 

N2LBF = 0.224809; % https://en.wikipedia.org/wiki/Newton_(unit)
IN2M = 0.0254;  % [m/in] Conversion factor from in to m
M2IN = 1 / IN2M;  % [in/m] Conversion factor from m to in

% Importing SFD and RFD
currentDirectory = pwd;
cd('C:../SFD')
sfdData = load('sfd_outputs_max_q.mat');
rfdData = load('rfd_outputs_recovery.mat');
lengthLoadssfd = sfdData.length_along_rocket_linspace * M2IN; % Length along the rocket converted to inches
lengthLoadsrfd = rfdData.length_along_rocket_linspace * M2IN; % Length along the rocket converted inches

override_loads = false;

if override_loads == false
    axialLoadssfd = sfdData.axial_array * N2LBF; % Axial loads converted to pounds
    shearLoadssfd = sfdData.shear_array * N2LBF; % Shear loads converted to pounds
    momentLoadssfd = sfdData.bending_array * 8.85; % Moment loads converted to inch-pounds
    
    axialLoadsrfd = rfdData.axial_array * N2LBF; % Axial loads converted to pounds
    shearLoadsrfd = rfdData.shear_array * N2LBF; % Shear loads converted to pounds
    momentLoadsrfd = rfdData.bending_array * 8.85; % Moment loads conevrted to inch-pounds
else
    axialLoadssfd = 200 * ones(1, numel(lengthLoadssfd));
    shearLoadssfd = 200 * ones(1, numel(lengthLoadssfd));
    momentLoadssfd = 30 * ones(1, numel(lengthLoadssfd));

    axialLoadsrfd = 200 * ones(1, numel(lengthLoadsrfd));
    shearLoadsrfd = 200 * ones(1, numel(lengthLoadsrfd));
    momentLoadsrfd = 30 * ones(1, numel(lengthLoadsrfd));
end

cd(currentDirectory)


%% Net Load Calcs

topLocation = max(location);
bottomLocation = min(location);

[~, locationMaxQ(1)] = min(abs(lengthLoadssfd - topLocation));
[~, locationMaxQ(2)] = min(abs(lengthLoadssfd - bottomLocation));
[~, locationRecovery(1)] = min(abs(lengthLoadsrfd - topLocation));
[~, locationRecovery(2)] = min(abs(lengthLoadsrfd - bottomLocation));

if radius == 0
    netCompressionMaxQ = axialLoadssfd; 
    netTensionMaxQ = axialLoadssfd; 
    netCompressionRecovery = axialLoadsrfd; 
    netTensionRecovery = axialLoadsrfd;
else
    netCompressionMaxQ = axialLoadssfd + (2 .* momentLoadssfd ./ radius); 
    netTensionMaxQ = axialLoadssfd - (2 .* momentLoadssfd / radius); 
    netCompressionRecovery = axialLoadsrfd + (2 .* momentLoadsrfd ./ radius); 
    netTensionRecovery = axialLoadsrfd - (2 .* momentLoadsrfd / radius);
end

maxMoment = max([momentLoadssfd(locationMaxQ(1)), momentLoadssfd(locationMaxQ(2)), momentLoadsrfd(locationRecovery(1)), momentLoadsrfd(locationRecovery(2))]);
maxAxial = max([axialLoadssfd(locationMaxQ(1)), axialLoadssfd(locationMaxQ(2)), axialLoadsrfd(locationRecovery(1)), axialLoadsrfd(locationRecovery(2))]);


max_Q_compressive_limit_load = max([netCompressionMaxQ(locationMaxQ(1)), netCompressionMaxQ(locationMaxQ(2))]);
recovery_compressive_limit_load = max([netCompressionRecovery(locationRecovery(1)), netCompressionRecovery(locationRecovery(2))]);

maxCompression = max([netCompressionMaxQ(locationMaxQ(1)), netCompressionRecovery(locationRecovery(1)), netCompressionMaxQ(locationMaxQ(2)), netCompressionRecovery(locationRecovery(2))]);
maxTension = min([netTensionMaxQ(locationMaxQ(1)), netTensionRecovery(locationRecovery(1)), netTensionMaxQ(locationMaxQ(2)), netTensionRecovery(locationRecovery(2))]);

%% GRAPHS

timeToGraph = exist('graphStatus');

if timeToGraph
    figure(1)
    subplot(1,3,1)
    plot(lengthLoadssfd, axialLoadssfd);
    hold on;
    plot([lengthLoadssfd(locationMaxQ(1)), lengthLoadssfd(locationMaxQ(2))], [axialLoadssfd(locationMaxQ(1)), axialLoadssfd(locationMaxQ(2))], '*r')
    hold off;
    grid on;
    title('Axial Loads at Max Q');
    xlabel('Length along Rocket (inches)');
    ylabel('Axial Load (pounds)');

    subplot(1,3,2)
    plot(lengthLoadssfd, shearLoadssfd);
    hold on;
    plot([lengthLoadssfd(locationMaxQ(1)), lengthLoadssfd(locationMaxQ(2))], [shearLoadssfd(locationMaxQ(1)), shearLoadssfd(locationMaxQ(2))], '*r')
    hold off;
    grid on;
    title('Shear Loads at Max Q');
    xlabel('Length along Rocket (inches)');
    ylabel('Shear Load (pounds)');

    subplot(1,3,3)
    plot(lengthLoadssfd, momentLoadssfd);
    hold on;
    plot([lengthLoadssfd(locationMaxQ(1)), lengthLoadssfd(locationMaxQ(2))], [momentLoadssfd(locationMaxQ(1)), momentLoadssfd(locationMaxQ(2))], '*r')
    hold off;
    grid on;
    title('Moment Loads at Max Q');
    xlabel('Length along Rocket (inches)');
    ylabel('Moment Load (pounds-in)');


    figure(2)
    subplot(1,3,1)
    plot(lengthLoadsrfd, axialLoadsrfd);
    hold on;
    plot([lengthLoadsrfd(locationRecovery(1)), lengthLoadsrfd(locationRecovery(2))], [axialLoadsrfd(locationMaxQ(1)), axialLoadsrfd(locationMaxQ(2))], '*r')
    hold off;
    grid on;
    title('Axial Loads during Recovery');
    xlabel('Length along Rocket (inches)');
    ylabel('Axial Load (pounds)');

    subplot(1,3,2)
    plot(lengthLoadsrfd, shearLoadsrfd);
    hold on;
    plot([lengthLoadsrfd(locationRecovery(1)), lengthLoadsrfd(locationRecovery(2))], [shearLoadsrfd(locationMaxQ(1)), shearLoadsrfd(locationMaxQ(2))], '*r')
    hold off;
    grid on;
    title('Shear Loads at Recovery');
    xlabel('Length along Rocket (inches)');
    ylabel('Shear Load (pounds)');

    subplot(1,3,3)
    plot(lengthLoadsrfd, momentLoadsrfd);
    hold on;
    plot([lengthLoadsrfd(locationRecovery(1)), lengthLoadsrfd(locationRecovery(2))], [momentLoadsrfd(locationMaxQ(1)), momentLoadsrfd(locationMaxQ(2))], '*r')
    hold off;
    grid on;
    title('Moment Loads at Recovery');
    xlabel('Length along Rocket (inches)');
    ylabel('Moment Load (pounds-in)');

end

%% For PDR Display Purposes

timeToPdr = exist('graphStatus');
timeToPdr = 1;
if timeToPdr
    for load_case = ["maxQ", "Recovery"]
        fprintf("Load Case: %s\n", load_case);
        for theta = [0, pi]
            fprintf("\tTheta: %.2f\n", theta);
            if (load_case == "maxQ")
                M = largest_value([momentLoadssfd(locationMaxQ(1)), momentLoadssfd(locationMaxQ(2))]);
                Fg = largest_value([axialLoadssfd(locationMaxQ(1)), axialLoadssfd(locationMaxQ(2))]);
            elseif (load_case == "Recovery")
                M = largest_value([momentLoadsrfd(locationRecovery(1)), momentLoadsrfd(locationRecovery(2))]);
                Fg = largest_value([axialLoadsrfd(locationRecovery(1)), axialLoadsrfd(locationRecovery(2))]);
            end

            R = radius;

            F3 = (Fg - ((2 * M * cos(theta)) / R)) / 3;
            F1 = M * cos(theta) / (3 * R) + Fg / 3 - M * sin(theta) / (R * sqrt(3));
            F2 = M * cos(theta) / (3 * R) + Fg / 3 + M * sin(theta) / (R * sqrt(3));
            fprintf("\t\tF1: %2.2f\n\t\tF2: %2.2f\n\t\tF3: %2.2f\n\n", F1, F2, F3)

        end
    end
end