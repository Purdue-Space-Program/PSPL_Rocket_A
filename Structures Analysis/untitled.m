%% AAE 204 HW3 Problem 3
%% David Gustafsson

applied_load_lbf = [ ...
    1000, 2000, 6000, 10000, 12000, 12900, 13400, 13600, ...
    13800, 14000, 14400, 15200, 16800, 18400, 20000, ...
    22400, 22600 ...
];

elongation_in = [ ...
    0.0002, 0.0006, 0.0019, 0.0033, 0.0039, 0.0043, ...
    0.0047, 0.0054, 0.0063, 0.0090, 0.0102, 0.0130, ...
    0.023, 0.0336, 0.0507, 0.1108, 0.12 ...
];

% Given initial values
initial_length_in   = 2;
initial_diameter_in = 0.505;

% Final minimum diameter
minimum_diameter_in = 0.42;

% Cross-sectional areas
initial_cross_sectional_area_in2 = pi * (initial_diameter_in / 2)^2;
minimum_cross_sectional_area_in2 = pi * (minimum_diameter_in / 2)^2;

% Engineering stress and strain
engineering_stress_psi = applied_load_lbf / initial_cross_sectional_area_in2;
engineering_strain     = elongation_in / initial_length_in;

% Plot stress-strain curve
plot(engineering_strain, engineering_stress_psi)
xlabel('Strain')
ylabel('Stress (psi)')
