import numpy as np
import matplotlib.pyplot as plt
import sfd
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import vehicle_parameters as vehicle
from scipy.io import savemat

location = "max_q" # Change to "max_q" or "off_the_rail"

LB2KG = 0.453592
FT2M = 0.3048
IN2M = 0.0254

burn_time = vehicle.parameters.burn_time # [s]
total_time = burn_time # [s]
cg_array, cg_max_q, cg_off_the_rail = sfd.updateCG(vehicle, burn_time, total_time) # Center of gravity over time

# plt.plot(cg_array)
# plt.show()

max_q_inputs = {"velocity": vehicle.parameters.max_velocity, "ax": vehicle.parameters.max_acceleration, "total_mass": vehicle.parameters.dry_mass}
off_the_rail_inputs = {"velocity": vehicle.parameters.off_the_rail_velocity, "ax": vehicle.parameters.off_the_rail_acceleration, "total_mass": vehicle.parameters.wet_mass}

location = "off_the_rail" # Change to "max_q" or "off_the_rail"

if location == "max_q":
    velocity = max_q_inputs["velocity"]
    ax = max_q_inputs["ax"]
    total_mass = max_q_inputs["total_mass"]
    mach = vehicle.parameters.max_mach
    cg = cg_max_q
elif location == "off_the_rail":
    velocity = off_the_rail_inputs["velocity"]
    ax = off_the_rail_inputs["ax"]
    total_mass = off_the_rail_inputs["total_mass"]
    mach = velocity / 343 # Speed of sound near sea level ~343 m/s
    cg = cg_off_the_rail

# Inputs
air_density = 1.81 # [kg / m^3] NEED
max_q_wind_gust = 13.4112 # [m / s] about 30 mph NEED
off_the_rail_rail_whip = 5 # [m / s] about 11 mph NEED
diameter = vehicle.parameters.tube_outer_diameter # [m]

if location == "max_q":
    linear_density_array, length_along_rocket_linspace = sfd.mass_model(vehicle.rocket_dict_dry)
    velocity = vehicle.parameters.max_velocity
    ax = vehicle.parameters.max_acceleration * 9.81
    mach = vehicle.parameters.max_mach
    cg = cg_max_q
    wind_gust = max_q_wind_gust
elif location == "off_the_rail":
    linear_density_array, length_along_rocket_linspace = sfd.mass_model(vehicle.rocket_dict_wet)
    velocity = vehicle.parameters.off_the_rail_velocity
    ax = vehicle.parameters.off_the_rail_acceleration * 9.81
    mach = velocity / 343 # Speed of sound near sea level ~343 m/s
    cg = cg_off_the_rail
    wind_gust = off_the_rail_rail_whip
total_length = length_along_rocket_linspace[-1] # [m]
total_mass = np.sum(linear_density_array * (length_along_rocket_linspace[1] - length_along_rocket_linspace[0])) # [kg]
print(f"Mass used: {total_mass} kg")
print(f"Length used: {total_length} m")

dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]  # [m]

thrust = vehicle.parameters.jet_thrust # [N]

# Fins
root_chord = 11 * IN2M # [m]
tip_chord = 2 * IN2M # [m]
sweep_length = (11 - 2) * IN2M # [m]
fin_height = 8 * IN2M # [m]
numFins = 3 # [m]
fin_top = vehicle.lower_fuel_bulkhead.bottom_distance_from_aft  # [m]
noseconeToFin = total_length - fin_top # [m]
print(f"noseconeToFin: {noseconeToFin} m")

# total_mass = 35.7 # [kg] This is when the bending moment will start and end at exactly 0 at max q
# total_mass = 35.3 # [kg] This is when the bending moment will start and end at exactly 0 at off the rail
# total_mass = np.sum(component.mass for component in vehicle.mass_distribution) # [kg] Total mass from mass distribution
# total_mass = np.sum(linear_density_array * (length_along_rocket_linspace[1] - length_along_rocket_linspace[0]))

# Calculated inputs
Q = sfd.calcQ(air_density, velocity)
AOA = sfd.calcAOA(wind_gust, velocity)
S = sfd.calcS(diameter)

# Calculated values
# Fins
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

# Converting to matlab file
matlab_dict = {"axial_array": axial_array, "shear_array": shear_array, "bending_array": bending_array, "length_along_rocket_linspace": length_along_rocket_linspace} # Dictionary to save as .mat file
if location == "max_q":
    savemat("sfd_outputs_max_q.mat", matlab_dict) # Save as .mat file for MATLAB
elif location == "off_the_rail":
    savemat("sfd_outputs_off_the_rail.mat", matlab_dict) # Save as .mat file for MATLAB

# Plotting
N2LBS = 0.224809
M2FT = 3.28084

plt.figure()
plot_num = 1
for variable in ["shear_array", "bending_array", "axial_array"]:
    if variable == "shear_array":
        plot = shear_array * N2LBS
        ylabel = "Shear Force [lbs]"
        title = f"Shear Forces at {location}"
    if variable == "bending_array":
        plot = bending_array * N2LBS
        ylabel = "Bending Moment [lbs-ft]"
        title = f"Bending Moments at {location}"
    if variable == "axial_array":
        plot = axial_array * N2LBS
        ylabel = "Axial Force [lbs]"
        title = f"Axial Forces at {location}"
    plt.subplot(1,3, plot_num)
    plt.plot(length_along_rocket_linspace * M2FT, plot)
    plt.title(title)
    plt.xlabel("Length from aft [ft]")
    plt.ylabel(ylabel)
    plt.grid()
    plot_num += 1
plt.show()

print(f"cg {cg:.2f} m")
print(f"finCP {finCP:.2f} m")
print(f"total_length {total_length:.2f} m")
print(f"Wet mass: {vehicle.parameters.wet_mass:.2f} kg")
print(f"Dry mass: {vehicle.parameters.dry_mass:.2f} kg")
print(f"Total mass using vehicle mass distribution: {np.sum(component.mass for component in vehicle.mass_distribution):.2f} kg")
print(f"Total mass using linear density array: {np.sum(linear_density_array * (length_along_rocket_linspace[1] - length_along_rocket_linspace[0])):.2f} kg")
print(vehicle.oxidizer_tank.bottom_distance_from_aft)
print(linear_density_array)
print(length_along_rocket_linspace)
print(f"ax: {ax} m/s^2")
print(f"rho: {air_density} kg/m^3")
print(f"V: {velocity} m/s")
print(f"Diameter: {S} m")
print(f"Angular acceleration: {r} rad/s^2")
'''
print("Parameters")
print(f"Air density: {air_density}")
print(f"Velocity: {velocity}")
print(f"Wind gust: {wind_gust}")
print(f"Diameter: {diameter}")

print(f"Root chord: {root_chord}")
print(f"Tip chord: {tip_chord}")
print(f"Sweep length: {sweep_length}")
print(f"Fin height: {fin_height}")
print(f"Number of fins: {numFins}")

print(f"Mach number: {mach}")

print(f"Burn time: {burn_time}")
print(f"Total time: {total_time}")

print(f"Nosecone to fin: {noseconeToFin}")
print(f"Total mass: {total_mass}")
print(f"Total length: {total_length}")
print(f"Thrust: {thrust}")
print(f"Max acceleration: {ax}")


print("Calculated values")
print(f"Dynamic pressure: {Q}")
print(f"Angle of attack: {AOA}")
print(f"Cross sectional area: {S}")

print(f"Fin stability derivative: {finSD}")
print(f"Mach coefficient: {machCoeff}")
print(f"Nose stability derivative: {noseSD}")
print(f"Nose lift: {noseLift}")
print(f"Fin lift: {finLift}")
print(f"Center of gravity: {cg}")
print(f"Inertia: {inertia}")
print(f"Lateral acceleration: {ay}")
print(f"Fin center of pressure: {finCP}")
print(f"Nose center of pressure: {noseCP}")
print(f"Angular acceleration: {r}")

'''