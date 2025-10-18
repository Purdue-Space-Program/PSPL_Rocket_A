import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sfd import getRocketSections, returnMassModel, getCG
from math import atan, pi

LB2KG = 0.453592
FT2M = 0.3048
IN2M = 0.0254

sfd_inputs = pd.ExcelFile('sfd_inputs.xlsx', engine='openpyxl')
rocket_dict = getRocketSections(sfd_inputs)
mass_model = returnMassModel(0.005)
totalMass = 0
for sec_name in rocket_dict:
    if 'mass' in rocket_dict[sec_name]:
        totalMass += rocket_dict[sec_name]['mass']
# print(totalMass) # TEST

# General parameters
diameter = 6 * IN2M # [m] Rocket cross-sectional diameter
S = pi * (diameter / 2)**2 # Rocket cross-sectional area

# Calculate Angle of attack
velWind = 1 # [m/s] Wind velocity # UPDATE
rocket_velocity = 1 # [m/s] Rocket velocity # UPDATE
AOA = atan(velWind / rocket_velocity) # [radians] Angle of attack

# Calculate dynamic pressure
rho = 1 # Fluid density # UPDATE
v = 1 # Fluid velocity # UPDATE
q = rho * v / 2

# Calculate lift-curve slopes
# dCLnose = 
# dCLfin = 
# dCLboat = 

# Calculate nosecone lift at AOA = 0
NNose = q * S * AOA * dCLnose


# Calculate nosecone lift at AOA = 0
def calculateNNose(Q, S, AOA, dCLnose):
    '''
    Q: Dynamic pressure
    S: Fuselage cross-sectional area
    AOA: Angle of attack
    dCLnose: Lift-curve slope of nose
    '''
    NNose = Q * S * AOA * dCLnose
    return NNose

# Calculate fin lift at AOA = 0
def calculateNFin(Q, S, AOA, dCLfin):
    NFin = Q * S * AOA * dCLfin
    return NFin

# Calculate boat tail lift at AOA = 0
def calculateLBoattail(Q, S, AOA, dCLboattail):
    NBoattail = (-1) * Q * S * AOA * dCLboattail
    return NBoattail

# Calculate angle of attack after encountering wind gust
def calculateAOA(velWind, velocity):
    AOA = atan(velWind / velocity)
    return AOA

# Calculate lateral acceleration after encountering wind gust
def calculateAy(NNose, NFin, NBoattail, totalMass):
    Ay = (NNose + NFin + NBoattail) / totalMass
    return Ay

# Calculate angular acceleration about vehicle center of gravity after encountering wind gust
def calculateR(NNose, NFin, NBoattail, MNose, MFin, MBoattail):
    '''
    MNose: Moment arm of nose from its center of pressure to CG
    MFin: Moment arm of fin from its center of pressure to CG
    MBoattail: Moment arm of boattail from its center of pressure to CG
    '''
    R = NNose * MNose + NFin * MFin + NBoattail * MBoattail
    return R

def calculateAxialLoad(T, DNose, DFuselage, DBase, S, deltaP):
    
    A = (-1) * T + DNose + DFuselage + DBase + S * deltaP