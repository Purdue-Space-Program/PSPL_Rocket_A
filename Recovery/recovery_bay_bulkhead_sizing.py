pi = 3.14159
import sys
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *
import vehicle_parameters as v

def shearStrength(material, diameter):
    return (material)*(pi*diameter*diameter*0.25)

def CalculateBearingStrength(material, diameter, thickness):
    bearing_strength = material * diameter * thickness
    return (bearing_strength)

def CalculateMoS(maximum_allowable_load, limit_load, FOS, fitting_factor):
    MoS = (maximum_allowable_load/(limit_load*FOS*fitting_factor)) - 1
    return (MoS)

# these are not necessary to be their own functions
# def shearUMOS(max, actual, uFOS, fitting_factor):
#     return (max/(actual*uFOS*fitting_factor)) - 1

# def bearingUMOS(max, actual, uFOS, fitting_factor):
#     return (max/(actual*uFOS*fitting_factor)) - 1

# def bearingYMOS(max, actual, sUMOS, bUMOS, uFOS, yFOS, fitting_factor):
#     if (sUMOS > bUMOS):
#         return (max/(actual*yFOS*fitting_factor)) - 1
#     else:
#         return (max/(actual*uFOS*fitting_factor)) - 1

def main():
    print("Recovery Bay Connector Bolted Joint")
    
    fSU_SS_316_U = 66000.0
    #fSU_AS_U = 180000.0
    #boltDiameter = 0.2075
    boltDiameter = 0.2614
    #boltDiameter = 0.32
    boltShearStrength = shearStrength(fSU_SS_316_U, boltDiameter)
    print(f"\tboltShearStrength: {boltShearStrength}")

    F_bru_A_6061_T6 = 88_000 # [psi] from MMPDS-2019: https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/1934065665/Common+Material+Properties
    F_bry_A_6061_T6 = 58_000 # [psi] from MMPDS-2019: https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/1934065665/Common+Material+Properties
    plateThickness = 0.125
    boltBearingStrengthUltimate = CalculateBearingStrength(F_bru_A_6061_T6, boltDiameter, plateThickness)
    boltBearingStrengthYield = CalculateBearingStrength(F_bry_A_6061_T6, boltDiameter, plateThickness)
    print(f"\tboltBearingStrengthUltimate: {boltBearingStrengthUltimate:.2f}")
    print(f"\tboltBearingStrengthYield: {boltBearingStrengthYield:.2f}")


    ultimateFOS = 2.0
    yieldFOS = 1.5
    
    #Need to update axial force.
    #netAxialForce = 300
    netAxialForce = 16230
    #numberOfBolts = 12
    numberOfBolts = 16
    shearForcePerBolt = netAxialForce/numberOfBolts

    if (boltShearStrength > boltBearingStrengthUltimate):
        print("\tBearing Critical!")
        fittingFOS = 1.15
    else: 
        print("\tSheer Critical!")
        fittingFOS = 2.0

    shearUltimateMOS = CalculateMoS(boltShearStrength, shearForcePerBolt, ultimateFOS, fittingFOS)
    bearingUltimateMOS = CalculateMoS(boltBearingStrengthUltimate, shearForcePerBolt, ultimateFOS, fittingFOS)
    bearingYieldMOS = CalculateMoS(boltBearingStrengthYield, shearForcePerBolt, yieldFOS, fittingFOS)
    print(f"\tshearUltimateMOS: {shearUltimateMOS:.2f}")
    print(f"\tbearingUltimateMOS: {bearingUltimateMOS:.2f}")
    print(f"\tbearingYieldMOS: {bearingYieldMOS:.2f}")

main()