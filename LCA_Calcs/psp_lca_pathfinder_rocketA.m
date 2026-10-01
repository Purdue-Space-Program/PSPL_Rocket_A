%--------------------------------------------------------------------------
% Chris Park 
% Lumped Capacitance Analysis
% Pathfinder Rocket A
%
% Calculates temperature change using lumped capacitance model.
% Requires material properties, dimensions, combustion temps, and 
% convection coefficients.
%--------------------------------------------------------------------------
clear
clc
close all

outputFolder = 'C:\Users\Chris\OneDrive\Desktop\PSP\MatLb\LCA Calcs\Outputs';
dataFolder = 'C:\Users\Chris\OneDrive\Desktop\PSP\MatLb\LCA Calcs\Inputs';

% IMPORTS
% Engine Contour
contour = readmatrix(fullfile(dataFolder, "chamber_contour_meters.csv"));
% Heat Transfer Coefficient
h_list = readmatrix(fullfile(dataFolder, "chamber_heat_transfer_coefficient.csv"));
% Gas Temp Along Wall
T_flame = readmatrix(fullfile(dataFolder, "chamber_ambient_temp.csv"));
% Pressure
wall_pressure = readmatrix(fullfile(dataFolder, "chamber_pressure.csv"));

% MATERIAL PROPERTIES
% 1018 Steel
steel_rho = 7870;           % [kg/m^3]      Density
steel_Cp  = 470;            % [J/(kg*K)]    Specific Heat
steel_k   = 51.9;           % [W/(m*K)]     Thermal Conductivity
steel_s   = 2.4e8;          % [Pa]          Yield Strength
steel_E   = 1.9e11;         % [Pa]          Young's Modulus
steel_v   = 0.29;           % [unitless]    Poisson's Ratio
steel_a   = 1.2e-5;         % [1/K]         Thermal Expansion

% BURN PARAMETERS
t_burn  = 2;                % [sec]         Burn Time
T_init  = 293;              % [K]           Initial Temperature

% SAFETY PARAMETERS
T_max         = 773;       % [K]            Maximum Safe Temperature
safety_factor = 1.25;

% --- CALCULATIONS ---

% Chamber Dimensions
x = contour(:,1);      % [m] Depth
r = contour(:,2);      % [m] Radius
N = numel(x);          % Number of data points
M = N - 1;             % Number of sections

% Surface Area
ds    = sqrt(diff(x).^2 + diff(r).^2);  % [m] Length of each segment
r_mid = (r(1:end-1) + r(2:end))/2;      % [m] Mean radius of each segment 
x_mid = (x(1:end-1) + x(2:end))/2;      % [m] Midpoint of each segment 
A_sec = 2*pi*r_mid .* ds;               % [m^2] Surface area

% Sectioning Values
h_sec = (h_list(1:end-1)+h_list(2:end))/2;          % [W/m^2-K] Sectioned h value
T_gas_wall = (T_flame(1:end-1)+T_flame(2:end))/2;   % [K] Section wall gas value

% - LCA Calculations -
log_term   = log((T_gas_wall - T_init) ./ (T_gas_wall - T_max));

Lc_min     = h_sec * t_burn ./ (steel_rho * steel_Cp .* log_term);  % [m] Minimum Required Thickness
Lc_build   = Lc_min * safety_factor;                                % [m] Design Thickness

% All further calculations done with the design thickness.

time_const = steel_rho * steel_Cp * Lc_build ./ h_sec;
T_final    = T_gas_wall + (T_init - T_gas_wall) .* exp(-t_burn ./ time_const); 

conduct = [];

for i = 1:length(T_final)
    conduct(i, 1)   = linearExtrap(T_final(i), steel_k);
end

% Biot Number for Validity
biot = h_sec .* Lc_build ./ conduct;
biot_index = min(find(biot > 0.1));

% [kg] Engine Mass 
sec_vol  = A_sec .* Lc_build;
sec_mass = steel_rho * sec_vol; 

% -- HOOP STRESS --
% Inputs
r_i        = r_mid;                  % [m]   Inner radius
t_wall     = Lc_build;               % [m]   Thickness
temp_dif   = T_gas_wall - T_final;   % [k]   Inner/Outer Wall Temp Diff
max_stress = [];                     % [Pa]  Maximum Allowable Stress
young_mod  = [];                     % [Pa]  Young's Modulus
expansion  = [];                     % [1/K] Thermal Expansion
poisson    = [];                     % []    Poisson's Ratio

for i = 1:length(T_final)
    young_mod(i, 1) = linearExtrap(T_final(i), steel_E);
    expansion(i, 1) = linearExtrap(T_final(i), steel_a);
    poisson(i, 1)   = linearExtrap(T_final(i), steel_v);
end

q_therm = h_sec .* temp_dif;
hoop_thermal = young_mod .* expansion .* q_therm .* t_wall ./ (2 * (0.71) .* conduct);

% Pressure Calculations
pressure      = (wall_pressure(1:end-1)+wall_pressure(2:end))/2;
hoop_pressure = pressure .* r_i ./ t_wall;

% Validity
hoop_valid = r_i ./ t_wall;
hoop_index = min(find(hoop_valid < 10));

% Stress Calculations
for u = 1:length(x_mid)
    max_stress(u, 1) = linearExtrap(T_final(u), steel_s);
end

rec_stress = max_stress ./ safety_factor; % [MPa] Target stress
tot_stress = hoop_thermal + hoop_pressure;

% --- PLOTTING ---

% Chamber Wall Thickness
figure(1); clf;
grid on; hold on;

plot(x(1:end-1), Lc_build, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'red');
plot(x(1:end-1), Lc_min, 'LineWidth', 1.5, 'Color', '#d41e11')
yline(max(Lc_build), 'Label', 'Maximum Thickness', 'Color', 'black','LineStyle','--')

xlabel('Depth [m]')
ylabel('Thickness [m]');
title('Chamber Wall Thickness');
legend('Recommended FoS: 1.25', 'Absolute Minumum', 'Location','southeast')

% Biot and Chamber Contour
figure(2); clf;
grid on; hold on;

plot(x_mid, biot, 'LineWidth', 2)
xline(x_mid(biot_index), 'LineWidth', 1.5, 'Label', 'Biot > 0.1', 'Color', '#4b0acc')
yline(0.1, '--r', "Max Valid Biot");

xlabel('Chamber Depth [m]')
ylabel('Biot Number [unitless]');
title('Biot Number Along Contour');

% Built Chamber Section
figure(3); clf;
hold on; grid on;

plot(x, r, 'Color', 'red', 'LineStyle','--')
plot(x_mid, r_mid + Lc_build, 'Color', 'black', 'LineWidth', 2)
plot(x, -r, 'Color', 'red', 'LineStyle','--')
plot(x_mid, -r_mid - Lc_build, 'Color', 'black', 'LineWidth', 2)
axis equal;

xlabel('Chamber Depth [m]')
ylabel('Chamber Radius [m]')
title('Chamber Cross Section')
legend('Inner Contour', 'Outer Contour')

% Hoop Stress
figure(4); clf;
hold on; grid on;

plot(x_mid, hoop_pressure, 'Color', '[0 0.2 0.9]', 'LineWidth', 1)
plot(x_mid, max_stress, 'Color', '[0 0.3 0.7]', 'LineWidth', 1, 'LineStyle','--')
plot(x_mid, rec_stress, 'Color', '[0 0.6 0.9]', 'LineWidth', 1, 'LineStyle','--')
xline(x_mid(hoop_index), 'Color', 'black', 'LineStyle', ':','Label','Thick to Thin Wall Hoop Stress')

title('Hoop Stress Over Chamber Contour')
xlabel('Chamber Depth [m]')
ylabel('Stress [Pa]')
legend('Hoop Stress', 'Max Stress', 'Target Stress', 'Location','southwest')

% Total Stress
figure(5); clf;
hold on; grid on;

plot(x_mid(1:end-150), tot_stress(1:end-150), 'Color', '[0.1 0.4 0.5]', 'LineWidth', 1)
plot(x_mid(1:end-150), max_stress(1:end-150), 'Color', '[0 0.3 0.7]', 'LineWidth', 1, 'LineStyle','--')
plot(x_mid(1:end-150), rec_stress(1:end-150), 'Color', '[0 0.6 0.9]', 'LineWidth', 1, 'LineStyle','--')
xline(x_mid(hoop_index),'Color', 'black', 'LineStyle', ':','Label','Thick to Thin Wall Hoop Stress')

title('Circumferential Stress Over Chamber Contour')
xlabel('Chamber Depth [m]')
ylabel('Stress [Pa]')
legend('Total Stress', 'Max Stress', 'Target Stress', 'Location','southwest')


% --- PRINTOUTS ---
fprintf("Max required thickness: %.2f mm\n", 1000*max(Lc_min));
fprintf("Max recommended thickness: %.2f mm\n\n", 1000*max(Lc_build));

fprintf("Max mean wall temp at end of burn with design thickness: %.0f K (limit %g K)\n", max(T_final), T_max);
fprintf("Total wall mass at design thickness: %.2f kg\n", sum(sec_mass));


% --- EXPORTS ---
exportgraphics(figure(5), fullfile(outputFolder, 'total_stress.png'))
exportgraphics(figure(4), fullfile(outputFolder, 'hoop_stress.png'))
exportgraphics(figure(3), fullfile(outputFolder, 'chamber_cross_section.png'))
exportgraphics(figure(2), fullfile(outputFolder, 'biot_number.png'))
exportgraphics(figure(1), fullfile(outputFolder, 'chamber_wall_thickness.png'))

% --- FUNCTIONS ---
function estim_val = linearExtrap(temp, input)
    % Physical Constant Temp Change
    dataFolder = 'C:\Users\Chris\OneDrive\Desktop\PSP\MatLb\LCA Calcs\Inputs';
    change_factor = readmatrix(fullfile(dataFolder, "physical_change.csv"));
    
    temps_range   = change_factor (:,1);
    related_val = change_factor (:,2);
    
    if temp < 200 || temp > 1200
        fprintf('Out of range')
        estim_val = -1;
        return
    else
        if temp >= 299 && temp < 700
            i = 1;
        else
            i = 2;
        end
    
        t0 = temps_range(i);      
        t1 = temps_range(i+1);
        y0 = related_val(i);   
        y1 = related_val(i+1);

        percent_factor = y0 + (y1 - y0) * (temp - t0) / (t1 - t0);
        
        estim_val = percent_factor * input;
    end
end