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

# print(inertia)
# print(r)
# print(mass_model)
shear_array = sfd.getShearForce(mass_model, totalLength, lift_dict, cp_dict, ay, r, cg)
bending_array = sfd.getBendingForce(shear_array, totalLength)
sfd.graphShear(shear_array, totalLength)
# graphBending(bending_array, totalLength)