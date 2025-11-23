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

air_density = 1.225 # [kg / m^3] NEED
velocity = vehicle.parameters.max_velocity # [m / s]
wind_gust = 9 # [m / s] about 69 mph NEED
diameter = vehicle.parameters.tube_outer_diameter # [m]

# Fins
root_chord = 0.29 # [m] NEED
tip_chord = 0.066 # [m] NEED
sweep_length = 0.224 # [m] NEED
fin_height = 0.136 # [m] NEED
numFins = 3 # [m] NEED

mach = vehicle.parameters.max_mach # [Mach number]

linear_density_array = vehicle.linear_density_array # [kg / m]
length_along_rocket_linspace = vehicle.length_along_rocket_linspace # [m]

burn_time = 2.24 # [s]
total_time = 5 # [s]

noseconeToFin = 2 # [m] NEED

total_mass = np.sum(component.mass for component in vehicle.mass_distribution) # [kg]
total_length = length_along_rocket_linspace[-1] # [m]

thrust = vehicle.parameters.jet_thrust # [N]
ax = vehicle.parameters.max_acceleration * 9.81 # [m / s]


Q = sfd.calcQ(air_density, velocity, wind_gust)
AOA = sfd.calcAOA(wind_gust, velocity)
S = sfd.calcS(diameter)

# Fins
finSD = sfd.calcFinSD(root_chord, tip_chord, sweep_length, fin_height, numFins, diameter)
machCoeff = sfd.calcMachCoeff(1, mach)
noseSD = sfd.calcNoseSD(machCoeff)
noseLift = sfd.calcLift(Q, S, AOA, noseSD)
finLift = sfd.calcLift(Q, S, AOA, finSD)
cg = sfd.calcCG(linear_density_array, length_along_rocket_linspace)
cg_array = sfd.updateCG(vehicle, burn_time, total_time)
inertia = sfd.calcRotationalInertia(linear_density_array, length_along_rocket_linspace, cg)
ay = sfd.calcLateralAcceleration(noseLift, finLift, total_mass)
finCP = sfd.calcFinCP(root_chord, tip_chord, sweep_length, fin_height, total_length, noseconeToFin)
noseCP = sfd.calcNoseCP(vehicle.nosecone.length, total_length)
r = sfd.calcAngularAcceleration(noseLift, finLift, noseCP, finCP, inertia, cg)
shear_array = np.array(sfd.calcShear(noseLift, finLift, noseCP, finCP, ay, linear_density_array, length_along_rocket_linspace, r, cg))
bending_array = np.array(sfd.calcBending(shear_array, length_along_rocket_linspace))
axial_array = np.array(sfd.calcAxial(thrust, ax, linear_density_array, length_along_rocket_linspace, air_density, 0.65, S, velocity)) # For medium size fins, Cd ~ 0.65 (UW Madison)
matlab_dict = {"axial_array": axial_array, "shear_array": shear_array, "bending_array": bending_array, "length_along_rocket_linspace": length_along_rocket_linspace}
savemat("sfd_outputs.mat", matlab_dict)

variable = "shear_array"
N2LBS = 0.224809
M2FT = 3.28084

if variable == "shear_array":
    plot = shear_array * N2LBS
    ylabel = "Shear Force [lbs]"
    title = "Shear Forces"
if variable == "bending_array":
    plot = bending_array * N2LBS
    ylabel = "Bending Moment [lbs-ft]"
    title = "Bending Moments"
if variable == "axial_array":
    plot = axial_array * N2LBS
    ylabel = "Axial Force [lbs]"
    title = "Axial Forces"
plt.plot(length_along_rocket_linspace * M2FT, plot)
plt.title(title)
plt.xlabel("Length from aft [ft]")
plt.ylabel(ylabel)
plt.grid()
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