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
point_masses = sfd.getPointMasses(rocket_dict)
mass_model = sfd.getMassModel(rocket_dict, point_masses, dy)
totalMass = sfd.getTotalMass(rocket_dict) # [kg]
totalLength = sfd.getTotalLength(rocket_dict) # [m]
cg = sfd.getCG(rocket_dict, point_masses, totalMass)

AOA = sfd.getAOA(aero_dict)
Q = sfd.getQ(aero_dict)
S = sfd.getArea(diameter) # Cross sectional area [m^2]

# Stability derivative
machCoeff = sfd.getMachAdjustedCoeff(1, 0.45)
noseSD = sfd.getNoseSD(machCoeff)
finSD = sfd.getFinSD(rocket_dict, diameter)

# Lift forces
noseLift = sfd.getLiftForce(Q, S, AOA, noseSD)
finLift = sfd.getLiftForce(Q, S, AOA, finSD)
boattailLift = 0
lift_dict = {'nose': noseLift, 'finLift': finLift, 'boattail': boattailLift}

# CP Location
noseconeToFin = totalLength - rocket_dict['fins']['length']
noseCP = sfd.getNoseCP(rocket_dict['nosecone']['length'], totalLength)
finCP = sfd.getFinCP(noseconeToFin, rocket_dict, totalLength)
boattailCP = 0
cp_dict = {'nose': noseCP, 'fin': finCP, 'boattail': boattailCP}

ay = sfd.getLatAccel(lift_dict, totalMass)
r = sfd.getAngularAccel(lift_dict, cp_dict, cg)

