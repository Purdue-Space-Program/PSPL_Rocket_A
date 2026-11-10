%% CHAMBER ANALYSIS AND WALL THICKNESS 

%% --- SETUP ---
clear
clc
close all

outputFolder = 'C:\Users\Chris\OneDrive\Desktop\PSP\MatLb\LCA Calcs\Outputs';
dataFolder = 'C:\Users\Chris\OneDrive\Desktop\PSP\MatLb\LCA Calcs\Inputs';

%% --- IMPORTS ---
% Engine Contour
contour = readmatrix(fullfile(dataFolder, "chamber_contour_meters.csv"));
% Heat Transfer Coefficient
h_list = readmatrix(fullfile(dataFolder, "chamber_heat_transfer_coefficient.csv"));
% Gas Temp Along Wall
T_flame = readmatrix(fullfile(dataFolder, "chamber_ambient_temp.csv"));
% Pressure
wall_pressure = readmatrix(fullfile(dataFolder, "chamber_pressure.csv"));
% Inner Wall Temp
wall_temp = readmatrix(fullfile(dataFolder, "chamber_temp.csv"));
% Heat Flux
q_list = readmatrix(fullfile(dataFolder, 'chamber_heat_flux.csv'));

%% --- CONSTANTS ---
% MATERIAL PROPERTIES
% 1018 Hot Rolled Steel
steel_rho = 7870;           % [kg/m^3]      Density
steel_Cp  = 470;            % [J/(kg*K)]    Specific Heat
steel_k   = 51.9;           % [W/(m*K)]     Thermal Conductivity
steel_s   = 2.4e8;          % [Pa]          Yield Strength
steel_E   = 1.9e11;         % [Pa]          Young's Modulus
steel_v   = 0.29;           % []            Poisson's Ratio
steel_a   = 1.2e-5;         % [1/K]         Thermal Expansion
steel_Ad  = 1.4e-5;         % [m^2/s]       Thermal Diffusivity

% BURN PARAMETERS
t_burn  = 2;                % [sec]         Burn Time
T_init  = 293;              % [K]           Initial Temperature

% SAFETY PARAMETERS
T_vec         = 900:-5:600; % [K]           Maximum Safe Temperature
elong_max     = 36;         % [%]           Max Elongation Percent
stress_max    = 2.2e8;      % [Pa]          Maximum Allowable Stress
safety_factor = 1.4;        % []            Are We Safe????

% OTHER NUMBERS
p_atm      = 101325;        % [Pa]          Atmospheric Pressure
pa_to_psi  = 0.000145038;   % [psi/Pa]      Conversion for Graphs

% DEFINED THINGS (random bullshit)
conduct    = [];
hoop_gas   = [];
hoop_therm = [];
biot       = [];
safe       = zeros(1, length(T_vec))';
count      = 1;

%% --- CALCULATIONS ---
% WALL THICKNESS
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
h_sec      = (h_list(1:end-1)+h_list(2:end))/2;     % [W/m^2-K] Sectioned h value
wall_sec   = (wall_temp(1:end-1)+wall_temp(2:end))/2;
T_gas_wall = (T_flame(1:end-1)+T_flame(2:end))/2;   % [K] Section wall gas value
pres_sec   = (wall_pressure(1:end-1)+wall_pressure(2:end))/2;
q_sec      = (q_list(1:end-1)+q_list(2:end))/2;

for T_max = T_vec
    [change_idx, biot, T_final, t_min, conduct] = lumpedCapacitance(h_sec, T_gas_wall, T_init, T_max, t_burn, steel_rho, steel_Cp, steel_k, safety_factor);
    

    r_o = t_min + r_mid;
    r_i = r_mid;

    b = r_o;
    a = r_i;

    % Pressure Hoop Stress
    gas_stress = pres_sec.*((r_i.^2+r_o.^2)./(r_o.^2-r_i.^2));

    % Temp Difference between inner and outer wall.
    dT = q_sec .* t_min ./ conduct;
    
    % Thermal Hoop Stress
    part_one = (-1 * linearExtrap(T_max, steel_a) * linearExtrap(T_max, steel_E) * dT) ./ (2*(1 - steel_v));
    part_two = (log(r_o./r_i) - 1)./log(r_o./r_i);
    part_three = (r_o.^2./r_i.^2 + 1)./(r_o.^2./r_i.^2 - 1);
    therm_stress = part_one .* (part_two + part_three);

    stress = abs(gas_stress + therm_stress);
    
    for i = 1:numel(stress)
        if (stress(i) < linearExtrap(T_max, stress_max))
            safe(count) = safe(count) + 1;
        end
    end

    biot(:,count) = h_sec .* t_min ./ conduct;

    count = count + 1;

    disp(T_max)
end

% [kg] Engine Mass
valid_num = 0;
for i = 1:length(safe)
    if (safe(i) > valid_num)
        valid_idx = i;
        valid_num = safe(i);
    end
end

[change_idx, biot, T_final, t_min, conduct] = lumpedCapacitance(h_sec, T_gas_wall, T_init, T_vec(valid_idx), t_burn, steel_rho, steel_Cp, steel_k, safety_factor);

sec_vol  = A_sec .* t_min * safety_factor;
sec_mass = steel_rho * sec_vol;

t_req = t_min * safety_factor;

%% --- OUTPUTS ---
fprintf("Max required thickness: %.2f mm\n", 1000*max(t_min));
fprintf("Max recommended thickness: %.2f mm\n\n", 1000*max(t_req));

fprintf("Max mean wall temp at end of burn with design thickness: %.0f K (limit %g K)\n", max(T_final), T_max);
fprintf("Total wall mass at design thickness: %.2f kg\n", sum(sec_mass));

%% --- EXPORTS ---

%% --- FUNCTIONS ---
function [biot_index, biot, T_final, Lc_min, conduct] = lumpedCapacitance(h_sec, T_gas_wall, T_init, T_max, t_burn, steel_rho, steel_Cp, steel_k, safety_factor)
    conduct = [];
    
    % LCA Calculations
    log_term   = log((T_gas_wall - T_init) ./ (T_gas_wall - T_max));
    
    Lc_min     = h_sec * t_burn ./ (steel_rho * steel_Cp .* log_term);  % [m] Minimum Required Thickness
        
    time_const = steel_rho * steel_Cp * Lc_min ./ h_sec;
    T_final    = T_gas_wall + (T_init - T_gas_wall) .* exp(-t_burn ./ time_const); 
    
    % Filling out conductivity
    for i = 1:length(T_final)
        conduct(i, 1)   = steel_k - 0.0333*(T_final(i)-273.15);
    end
    
    % Biot Number for Validity
    biot = h_sec .* Lc_min ./ conduct;
    biot_index = min(find(biot > 0.1));
end


function estim_val = linearExtrap(temp, input)
    % Physical Constant Temp Change
    dataFolder = 'C:\Users\Chris\OneDrive\Desktop\PSP\MatLb\LCA Calcs\Inputs';
    change_factor = readmatrix(fullfile(dataFolder, "physical_change.csv"));
    
    temps_range   = change_factor (:,1);
    related_val = change_factor (:,2);
     
    if temp < 200 || temp > 1200
        fprintf('Out of range bruh\n')
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

%% --- PLOTS ---
% Chamber Wall Thickness
figure(1); clf;
grid on; hold on;

plot(x_mid, t_req, 'Color', 'black', 'LineWidth', 2)
plot(x_mid, t_min, 'LineWidth', 1.5, 'Color', '#d41e11')
xline(x_mid(change_idx), 'Label', 'Invalid Biot', 'Color', 'red', LineStyle='--')
yline(max(t_req), 'Label', 'Maximum Thickness', 'Color', 'black','LineStyle','--')

xlabel('Depth [m]')
ylabel('Thickness [m]');
title('Chamber Wall Thickness');
legend('Recommended FoS: 1.4', 'Absolute Minumum', 'Location','southeast')

% Built Chamber Section
figure(2); clf;
hold on; grid on;

plot(x, r, 'Color', 'red', 'LineStyle','--')
plot(x_mid, r_mid + t_req, 'Color', 'black', 'LineWidth', 2)
plot(x, -r, 'Color', 'red', 'LineStyle','--')
plot(x_mid, -r_mid - t_req, 'Color', 'black', 'LineWidth', 2)
axis equal;

xlabel('Chamber Depth [m]')
ylabel('Chamber Radius [m]')
title('Chamber Cross Section')
legend('Inner Contour', 'Outer Contour')

% Hoop Stress Only
figure(3); clf;
hold on; grid on;

plot(x_mid, gas_stress * pa_to_psi, 'b')
plot(x_mid, therm_stress * pa_to_psi, color=[0.5 0.1 0.2])
plot(x_mid, stress * pa_to_psi, color='magenta')
yline(linearExtrap(T_vec(valid_idx), stress_max) * pa_to_psi, 'Label', 'Yield Stress', LineStyle='--', Color='r')
title('Hoop Stress Over the Chamber')
xlabel('Chamber Depth[m]')
ylabel('Stress [psi]')
legend("Physical Stress", "Thermal Stress", "Total Stress", 'Location','northeast')

% Inputs
figure(4); clf;

subplot(2,2,1)
hold on; grid on;

plot(x, r, 'black')
plot(x, -r, 'black')

title('Chamber Contour')
xlabel('Chamber Depth [m]')
ylabel('Chamber Radius [m]')
axis equal;

subplot(2, 2, 2)
hold on; grid on;

plot(x, h_list, 'red')

title('Convective Heat Transfer Coefficient')
xlabel('Chamber Depth [m]')
ylabel('Convective Heat Transfer [W/(m^2*K)]')

subplot(2, 2, 3)
hold on; grid on;

plot(x, T_flame, 'Color', '#FFA500')

title('Inner Wall Temperature')
xlabel('Chamber Depth [m]')
ylabel('Gas Temperature [K]')

subplot(2, 2, 4)
hold on; grid on;

plot(x, wall_pressure * pa_to_psi, 'blue')

title('Wall Gas Pressure Along Contour')
xlabel('Chamber Depth [m]')
ylabel('Pressure [psi]')

figure(5); clf;
hold on; grid on;

plot(flip(T_vec), safe, 'rs', 'LineWidth', 2)
title('Safe Points as a Function of Max Temperature of LCA')
xlabel('Max Allowable Temperature [K]')
ylabel('Number of Safe Points')

exportgraphics(figure(5), fullfile(outputFolder, '2safe_points.png'))
exportgraphics(figure(4), fullfile(outputFolder, '2chamber_info.png'))
exportgraphics(figure(3), fullfile(outputFolder, '2total_stress.png'))
exportgraphics(figure(2), fullfile(outputFolder, '2built_section.png'))
exportgraphics(figure(1), fullfile(outputFolder, '2wall_thickness.png'))
