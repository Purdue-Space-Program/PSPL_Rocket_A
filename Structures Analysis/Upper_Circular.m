% Upper Strut Sizing Code (ASD-based per ADM 2020)
% Dylan Agnese
% === Inputs ===
length = 1.614*12; % length of upper strut (in)
E = 1.0298e7; % Young's Modulus of Aluminum 6061 (psi)
yield_strength = 40611; % Yield strength (psi)
G = 3.77e6; % Shear modulus (psi)
density = 0.0975; % Material density (lb/in^3)
axial_load = 1217; % Axial compression (lb)
tension_load = 1009.852; % Max tensile load (lb)
bending = 0; % Bending moment (in-lb)
radius_rocket = (8.625/2)*1/12 - (1.1/2)*(1/12); % Rocket radius (in)
strut_od = .75; % Outer diameter (in)
strut_id = .58; % Inner diameter (in)
K = 0.9; % Effective length factor (assumed)
v = 0.33; %poisson ratio

% === Geometry Calculations ===
thickness = (strut_od - strut_id) / 2; % Wall thickness
area_total = (pi/4) * (strut_od^2 - strut_id^2); % Cross-sectional area
I = (pi/64) * (strut_od^4 - strut_id^4); % Moment of inertia
radius_gyration = sqrt(I / area_total); % Radius of gyration
slenderness_ratio = (K * length) / radius_gyration; % Slenderness ratio
transitional_ratio = sqrt((2*pi^2*E)/(K^2*yield_strength));
wt_ratio = strut_od / thickness; % Width-to-thickness ratio
plastic_mod = (pi/4) * ((strut_od/2)^3 - (strut_id/2)^3); % Plastic modulus
section_modulus = (pi * (strut_od^4 - strut_id^4)) / (32 * strut_od); % Section modulus


% Buckling Constants%

if wt_ratio <= ((0.11 * E) / yield_strength)
    Area  = area_total;
    Fe = (pi^2 * E) / (slenderness_ratio ^2);
    fprintf("Nonslender – local buckling not critical.\n");
else 
    Area = area_total * ((0.038* E) / (yield_strength *(wt_ratio)) + 2/3);
    Fe = (pi^2 * E) / (slenderness_ratio ^2);
    fprintf("Slender – local buckling must be considered.\n");
end

if slenderness_ratio <= (4.71 * (E / yield_strength)^ (1/2))
    Fcr = ((0.658) ^ (yield_strength / Fe)) * yield_strength;
else
    Fcr = 0.877 * Fe;
end






% === Allowable Strengths (ASD) ===
Omega_c = 1.65; % Compression safety factor (ADM ASD)
Omega_b = 1.67; % Bending safety factor (ADM ASD)

P_allow = (Fcr * Area) / Omega_c; % Allowable axial load (lb)

% === Combined Axial and Bending Check ===
area_total = (pi/4) * (strut_od^2 - strut_id^2);
mass = area_total * length * density;



% === Tension Check ===
tensile_stress = tension_load / area_total;
FoS_tension = yield_strength / tensile_stress;
 
% === Output ===
fprintf("Available Axial Strength (ASD): %.2f lbs\n", P_allow);
fprintf("Mass of one strut: %.2f lbs\n", mass);
fprintf("Strut OD: %.2f in, ID: %.2f in\n", strut_od, strut_id);
fprintf("Tensile FoS: %.2f\n", FoS_tension);
