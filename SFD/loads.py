import numpy as np
import matplotlib.pyplot as plt
import sfd
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import vehicle_parameters as vehicle
import constants as c
from scipy.io import savemat
import parseWind as pw
# ------------------------------------------------------------------------------

# Constants
LB2KG = 0.453592 # Pounds to Kilograms
FT2M = 0.3048 # Feet to Meters
IN2M = 0.0254 # Inches to Meters
N2LBF = 0.224809 # Newtons to Pounds force
gravity = 9.81 # [m / s^2]
# ------------------------------------------------------------------------------

# Select location to analyze
location = "max_q" # Change to "max_q" or "off_the_rail"
plot_on = False # Set to True to plot results, False to not plot
# ------------------------------------------------------------------------------

# Inputs
diameter = vehicle.parameters.tube_outer_diameter # [m]
thrust = vehicle.parameters.jet_thrust # [N]
total_length = vehicle.parameters.total_length # [m]

air_density = 1.81 # [kg / m^3] NEED
max_q_wind_gust = pw.percentile_75_wind_speed # [m / s] about 30 mph NEED
off_the_rail_rail_whip = 5 # [m / s] about 11 mph NEED

if location == "max_q":
    velocity = vehicle.parameters.six_DoF_max_velocity # [m / s] Velocity at max q
    ax = vehicle.parameters.one_DoF_max_acceleration * gravity # [m / s^2] Acceleration at max q
    total_mass = vehicle.parameters.dry_mass # [kg] Mass at max q
    mach = vehicle.parameters.six_DoF_max_mach # [] Mach number at max q
    cg = vehicle.parameters.dry_COM_location_from_bottom # [m] Center of gravity from bottom of rocket
    linear_density_array, length_along_rocket_linspace = sfd.mass_model(vehicle.rocket_dict_dry)
    wind_gust = max_q_wind_gust # [m / s] Wind gust at max q
elif location == "off_the_rail":
    velocity = vehicle.parameters.six_DoF_off_the_rail_velocity # [m / s] Velocity off the rail
    ax = vehicle.parameters.six_DoF_off_the_rail_acceleration * gravity # [m / s^2] Acceleration off the rail
    total_mass = vehicle.parameters.wet_mass # [kg] Mass off the rail
    mach = velocity / 343 # Speed of sound near sea level ~343 m/s
    cg = vehicle.parameters.wet_COM_location_from_bottom # [m] Center of gravity from bottom of rocket
    linear_density_array, length_along_rocket_linspace = sfd.mass_model(vehicle.rocket_dict_wet)
    wind_gust = off_the_rail_rail_whip # [m / s] Rail whip off the rail
else:
    raise ValueError("Invalid location selected. Choose 'max_q' or 'off_the_rail'.")

total_mass = np.sum(linear_density_array * (length_along_rocket_linspace[1] - length_along_rocket_linspace[0])) # [kg]

dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]  # [m]

# Fins
root_chord = 11 * IN2M # [m]
tip_chord = 2 * IN2M # [m]
sweep_length = (11 - 2) * IN2M # [m]
fin_height = 8 * IN2M # [m]
numFins = 3 # [m]
fin_top = vehicle.lower_fuel_bulkhead.bottom_distance_from_aft # [m]
noseconeToFin = total_length - fin_top # [m]

# Calculated inputs
Q = sfd.calcQ(air_density, velocity)
AOA = sfd.calcAOA(wind_gust, velocity) # [radians] # 1.02331 * np.pi / 180 # NEED
S = sfd.calcS(diameter) # [m^2] Cross sectional area
# ------------------------------------------------------------------------------

# Calculated values
finCP = sfd.calcFinCP(root_chord, tip_chord, sweep_length, fin_height, total_length, noseconeToFin) # Fin center of pressure
noseCP = sfd.calcNoseCP(vehicle.nosecone.length, total_length) # Nose center of pressure
finSD = sfd.calcFinSD(root_chord, tip_chord, sweep_length, fin_height, numFins, diameter) # Fin stability derivative
machCoeff = sfd.calcMachCoeff(1, mach) # Mach coefficient
noseSD = sfd.calcNoseSD(cg, noseCP, diameter) # Nose stability derivative
noseLift = sfd.calcLift(Q, S, AOA, noseSD) # Nose lift
finLift = sfd.calcLift(Q, S, AOA, finSD) # Fin lift
inertia = sfd.calcRotationalInertia(linear_density_array, length_along_rocket_linspace, cg) # Rotational inertia
ay = sfd.calcLateralAcceleration(noseLift, finLift, total_mass) # Lateral acceleration
r = sfd.calcAngularAcceleration(noseLift, finLift, noseCP, finCP, inertia, cg) # Angular acceleration
shear_array = np.array(sfd.calcShear(noseLift, finLift, noseCP, finCP, ay, linear_density_array, length_along_rocket_linspace, r, cg)) # Shear force array
bending_array = np.array(sfd.calcBending(shear_array, length_along_rocket_linspace)) # Bending moment array
axial_array = np.array(sfd.calcAxial(thrust, ax, linear_density_array, length_along_rocket_linspace, air_density, 0.65, S, velocity)) # Axial forces array, For medium size fins, Cd ~ 0.65 (UW Madison)
# ------------------------------------------------------------------------------

# Converting to matlab file
matlab_dict = {f"axial_array": axial_array, f"shear_array": shear_array, f"bending_array": bending_array, "length_along_rocket_linspace": length_along_rocket_linspace} # Dictionary to save as .mat file
if location == "max_q":
    savemat("sfd_outputs_max_q.mat", matlab_dict) # Save as .mat file for MATLAB
elif location == "off_the_rail":
    savemat("sfd_outputs_off_the_rail.mat", matlab_dict) # Save as .mat file for MATLAB
# ------------------------------------------------------------------------------

# Outputs
print(f"Outputs at {location}:")
print(f"Max shear force at {location}: {max(shear_array) * N2LBF:.2f} lbf")
print(f"Max bending moment at {location}: {max(bending_array) * N2LBF * c.M2FT:.2f} lbf-ft")
print(f"Max axial force at {location}: {max(axial_array) * N2LBF:.2f} lbf")
print("-----------------------------------")

print("Inputs:")
print(f"Wind gust at {location}: {wind_gust} m/s")
print(f"Air density at {location}: {air_density} kg/m^3")
print(f"Angle of attack at {location}: {AOA * (180 / np.pi):.2f} degrees")
print(f"Velocity at {location}: {velocity} m/s")
print(f"Acceleration at {location}: {ax} m/s^2")
print(f"Rocket mass at {location}: {total_mass:.2f} kg")
print(f"Center of gravity at {location}: {cg:.2f} m from bottom")
print(f"Cross sectional area: {S:.4f} m^2")
print("-----------------------------------")

print("Fin parameters:")
print(f"Fin center of pressure: {finCP:.2f} m from bottom")
print(f"Fin stability derivative: {finSD:.4f} ")
print(f"Fin lift: {finLift:.2f} N")
print(f"Root chord: {root_chord:.2f} m")
print(f"Tip chord: {tip_chord:.2f} m")
print(f"Sweep length: {sweep_length:.2f} m")
print(f"Fin height: {fin_height:.2f} m")
print(f"Number of fins: {numFins}")
print("-----------------------------------")

print("Nosecone parameters:")
print(f"Nose center of pressure: {noseCP:.2f} m from bottom")
print(f"Nose stability derivative: {noseSD:.4f} ")
print(f"Nose lift: {noseLift:.2f} N")
print("-----------------------------------")

# Plotting
plt.figure()
plot_num = 1
for variable in ["shear_array", "bending_array", "axial_array"]:
    if variable == "shear_array":
        plot = shear_array * N2LBF
        ylabel = "Shear Force [lbf]"
        title = f"Shear Forces at {location}"
    if variable == "bending_array":
        plot = bending_array * N2LBF * c.M2FT
        ylabel = "Bending Moment [lbf-ft]"
        title = f"Bending Moments at {location}"
    if variable == "axial_array":
        plot = axial_array * N2LBF
        ylabel = "Axial Force [lbf]"
        title = f"Axial Forces at {location}"
    plt.subplot(1,3, plot_num)
    plt.plot(length_along_rocket_linspace * c.M2FT, plot)
    plt.title(title)
    plt.xlabel("Length from aft [ft]")
    plt.ylabel(ylabel)
    plt.grid()
    plot_num += 1
if plot_on:
    plt.show()

