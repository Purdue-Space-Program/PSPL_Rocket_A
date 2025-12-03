import numpy as np
import matplotlib.pyplot as plt
import sfd
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import vehicle_parameters as vehicle
from scipy.io import savemat


LB2KG = 0.453592
FT2M = 0.3048
IN2M = 0.0254

burn_time = vehicle.parameters.burn_time # [s]
total_time = burn_time # [s]
cg_array, cg_max_q, cg_off_the_rail = sfd.updateCG(vehicle, burn_time, total_time) # Center of gravity over time

max_q_inputs = {"velocity": vehicle.parameters.max_velocity, "ax": vehicle.parameters.max_acceleration, "total_mass": vehicle.parameters.dry_mass}
off_the_rail_inputs = {"velocity": vehicle.parameters.off_the_rail_velocity, "ax": vehicle.parameters.off_the_rail_acceleration, "total_mass": vehicle.parameters.wet_mass}

location = "max_q" # Change to "max_q" or "off_the_rail"

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
air_density = 1.225 # [kg / m^3] NEED
max_q_wind_gust = 9 # [m / s] about 69 mph NEED
off_the_rail_wind_gust = 5 # [m / s] about 45 mph NEED
diameter = vehicle.parameters.tube_outer_diameter # [m]

# Fins
root_chord = 0.29 # [m] NEED
tip_chord = 0.066 # [m] NEED
sweep_length = 0.224 # [m] NEED
fin_height = 0.136 # [m] NEED
numFins = 3 # [m] NEED

linear_density_array = vehicle.linear_density_array # [kg / m]
length_along_rocket_linspace = vehicle.length_along_rocket_linspace # [m]


noseconeToFin = 2 # [m] NEED

total_mass = np.sum(component.mass for component in vehicle.mass_distribution) # [kg]
total_length = length_along_rocket_linspace[-1] # [m]

thrust = vehicle.parameters.jet_thrust # [N]
ax = vehicle.parameters.max_acceleration * 9.81 # [m / s]


Q = sfd.calcQ(air_density, velocity, max_q_wind_gust)
AOA = sfd.calcAOA(max_q_wind_gust, velocity)
S = sfd.calcS(diameter)

# Calculated values
# Fins
finSD = sfd.calcFinSD(root_chord, tip_chord, sweep_length, fin_height, numFins, diameter) # Fin stability derivative
machCoeff = sfd.calcMachCoeff(1, mach) # Mach coefficient
noseSD = sfd.calcNoseSD(machCoeff) # Nose stability derivative
noseLift = sfd.calcLift(Q, S, AOA, noseSD) # Nose lift
finLift = sfd.calcLift(Q, S, AOA, finSD) # Fin lift
cg = sfd.calcCG(linear_density_array, length_along_rocket_linspace) # Center of gravity
cg_array = sfd.updateCG(vehicle, burn_time, total_time) # Center of gravity over time
cg = cg_array[-1] # Final center of gravity
inertia = sfd.calcRotationalInertia(linear_density_array, length_along_rocket_linspace, cg) # Rotational inertia
ay = sfd.calcLateralAcceleration(noseLift, finLift, total_mass) # Lateral acceleration
finCP = sfd.calcFinCP(root_chord, tip_chord, sweep_length, fin_height, total_length, noseconeToFin) # Fin center of pressure
noseCP = sfd.calcNoseCP(vehicle.nosecone.length, total_length) # Nose center of pressure
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