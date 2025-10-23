import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sfd

LB2KG = 0.453592
FT2M = 0.3048
IN2M = 0.0254

# Parameters
diameter = 6 * IN2M
dy = 0.005

sfd_inputs = pd.ExcelFile('sfd_inputs.xlsx', engine='openpyxl')
rocket_dict = sfd.getRocketSections(sfd_inputs)
aero_dict = sfd.getAeroProperties(sfd_inputs)
point_masses = sfd.getPointMasses(sfd_inputs)
mass_model = sfd.getMassModel(rocket_dict, point_masses, dy)
totalMass = sfd.getTotalMass(rocket_dict) # [kg]
totalLength = sfd.getTotalLength(rocket_dict) # [m]
cg = sfd.getCG(rocket_dict, point_masses, totalMass)

location = ['off_the_rail', 'max_q']
AOA = sfd.getAOA(aero_dict, location[1])
Q = sfd.getQ(aero_dict, location[1])
S = sfd.getArea(diameter) # Cross sectional area [m^2]

# Stability derivative
machCoeff = sfd.getMachAdjustedCoeff(1, 0.45)
noseSD = sfd.getNoseSD(machCoeff)
finSD = sfd.getFinSD(rocket_dict, diameter)

# Lift forces
noseLift = sfd.getLiftForce(Q, S, AOA, noseSD)
finLift = sfd.getLiftForce(Q, S, AOA, finSD)
boattailLift = 0
lift_dict = {'nose': noseLift, 'fin': finLift, 'boattail': boattailLift}

# CP Location
noseconeToFin = totalLength - rocket_dict['fins']['length'] # Distance from nosecone to fin
noseCP = sfd.getNoseCP(rocket_dict['nosecone']['length'], totalLength) # Nose center of pressure
finCP = sfd.getFinCP(noseconeToFin, rocket_dict, totalLength) # Fin center of pressure
boattailCP = 0 # Boattail center of pressure
cp_dict = {'nose': noseCP, 'fin': finCP, 'boattail': boattailCP} # Center of pressure dictionary

inertia = sfd.getRotationalInertia(mass_model, cg, totalLength) # Inertia
ay = sfd.getLatAccel(lift_dict, totalMass) # Lateral acceleration
r = sfd.getAngularAccel(lift_dict, cp_dict, cg, inertia) # Angular acceleration

def getShearForce(mass_model, totalLength, lift_dict, cp_dict, ay, r):
    dy = totalLength / len(mass_model)
    lengths = [(i + 1) * dy for i in range(len(mass_model))]
    print(lengths) # TEST
    
    shear_array = (-1) * (ay * np.cumsum(mass_model) + r * np.cumsum(mass_model * (cg - lengths)))
    shear_array[int(cp_dict['nose'] / dy):] += lift_dict['nose']
    shear_array[int(cp_dict['fin'] / dy):] += lift_dict['fin']
    shear_array[int(cp_dict['boattail'] / dy):] -= lift_dict['boattail']
    
    return shear_array

def getBendingForce(shear_array, totalLength):
    dy = totalLength / len(shear_array)
    bending_array = np.cumsum(shear_array) * dy
    return bending_array

def graphShear(shear_array, totalLength):
    dx = totalLength / len(shear_array)
    x = [i * dx for i in range(len(shear_array))]
    plt.plot(x, shear_array)
    plt.title("Shear Forces")
    plt.show()

def graphBending(bending_array, totalLength):
    dx = totalLength / len(bending_array)
    x = [i * dx for i in range(len(bending_array))]
    plt.plot(x, bending_array)
    plt.title("Bending Forces")
    plt.show()

# print(inertia)
# print(r)
getShearForce(mass_model, totalLength, lift_dict, cp_dict, ay, r)
print(totalLength)
# print(mass_model)
print(cp_dict['fin'])
shear_array = getShearForce(mass_model, totalLength, lift_dict, cp_dict, ay, r)
bending_array = getBendingForce(shear_array, totalLength)
# graphShear(shear_array, totalLength)
# graphBending(bending_array, totalLength)