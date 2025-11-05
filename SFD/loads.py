import numpy as np
import matplotlib.pyplot as plt
import sfd
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import vehicle_parameters as vehicle

LB2KG = 0.453592
FT2M = 0.3048
IN2M = 0.0254

air_density = 1.225 # [kg / m^3] NEED
velocity = vehicle.parameters.max_velocity # [m / s]
wind_gust = 30 # [m / s] about 69 mph NEED
diameter = vehicle.parameters.tube_outer_diameter # [m]

# Fins
root_chord = 0.3 # [m] NEED
tip_chord = 0.1 # [m] NEED
sweep_length = 0.2 # [m] NEED
fin_height = 0.3 # [m] NEED
numFins = 3 # [m] NEED

mach = vehicle.parameters.max_mach

linear_density_array = vehicle.linear_density_array # [kg / m]
length_along_rocket_linspace = vehicle.length_along_rocket_linspace # [m]

burn_time = 2.24 # [s]
total_time = 5 # [s]

noseconeToFin = 2 # [m] NEED

total_mass = np.sum(component.mass for component in vehicle.mass_distribution) # [kg]
total_length = length_along_rocket_linspace[-1] # [m]

thrust = vehicle.parameters.jet_thrust
ax = vehicle.parameters.max_acceleration


Q = sfd.calcQ(air_density, velocity)
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
axial_array = np.array(sfd.calcAxial(thrust, ax, linear_density_array, length_along_rocket_linspace))


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
'''
plt.plot(length_along_rocket_linspace * M2FT, plot)
plt.title(title)
plt.xlabel("Length from aft [ft]")
plt.ylabel(ylabel)
plt.show()

'''
print(finSD)
print(noseSD)