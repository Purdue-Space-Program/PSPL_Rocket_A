%% Contants

LiftCoef = 0.377; % Taken from RAS Aero


%Fin dimensions
SweepDist = 10; % [in]
RootChord = 16; % [in]
tipChord = 6; % [in] 
span = 7.6; % [in]

maxMach = 0.4; % Taken from TrajArray
maxMachAlt = 2343; % [ft] Taken from trajArray

%% Calculations
% Pressure and temperature calculations (NASA Earth Atmosphere Model Formulas)
[~, speedSound, ~ , density] = atmosisa(maxMachAlt * 0.3048);

% Unit conversions
%speedSound = speedSound; % [m/s]
%densityUSA = density; % [kg/m^3]

trueVelo = speedSound * maxMach; % [m/s]
%trueVelo = maxVelo;

finArea_in = ((0.5 * SweepDist * span) + (tipChord * span) + (0.5 * (RootChord - SweepDist - tipChord) * span)); %[in^2]
finArea_m = finArea_in/144*0.09290304;

lift = LiftCoef * 0.5 * finArea_m * trueVelo^2 * density; %[N] Total lift on the wing

wing_pres = lift/finArea_m; % pressure calculation
fprintf("Wing pressure: %0.2f N/m\n", wing_pres)

liftUSA = lift * 0.2248089431; % [lbf]

fprintf("Lift: %.2fN = %0.2flbf\n", lift, liftUSA);

xStep = 0.00001;

countersweepDist = RootChord - tipChord - SweepDist;
sweepAngle = atan(SweepDist/span);
countersweepAngle = atan(countersweepDist/span);

Torque = 0;

span = span*0.0254; % conversion to meters
RootChord = RootChord*0.0254; % conversion to meters

for x = 0:xStep:span
    dA = (RootChord - (tan(sweepAngle) * x) - (tan(countersweepAngle) * x)) * xStep;
    dTorque = x * dA * wing_pres; % determine torque for slice
    Torque = Torque + dTorque; % sum torque
end

Torque_in = Torque*8.8507458; % lbf-in

fprintf("Total Torque = %.2f lb-in\n", Torque_in)

Torque_ft = Torque*0.737562; % lbf-ft

fprintf("Total Torque = %.2f lb-ft\n", Torque_ft)
