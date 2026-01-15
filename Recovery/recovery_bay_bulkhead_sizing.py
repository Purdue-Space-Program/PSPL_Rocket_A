pi = 3.14159
import sys
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *
import vehicle_parameters as v

def shearStrength(material, diameter):
    return (material)*(pi*diameter*diameter*0.25)

def bearingStrength(material, diameter, thickness):
    return (material)*(diameter)*(thickness)

def shearUMOS(max, actual, uFOS, fFOS):
    return (max/(actual*uFOS*fFOS)) - 1

def bearingUMOS(max, actual, uFOS, fFOS):
    return (max/(actual*uFOS*fFOS)) - 1

def bearingYMOS(max, actual, yFOS, fFOS):
    return (max/(actual*yFOS*fFOS)) - 1

def main():
    fSU_SS_316_U = 66000.0
    #boltDiameter = 0.2075
    boltDiameter = 0.32
    boltShearStrength = shearStrength(fSU_SS_316_U, boltDiameter)
    print(f"boltShearStrength: {boltShearStrength}")

    fBRU_A_6061_U = 88000.0
    fBRU_A_6061_Y = 58000.0
    plateThickness = 0.125
    boltBearingStrengthUltimate = bearingStrength(fBRU_A_6061_U, boltDiameter, plateThickness)
    boltBearingStrengthYield = bearingStrength(fBRU_A_6061_Y, boltDiameter, plateThickness)
    print(f"boltBearingStrengthUltimate: {boltBearingStrengthUltimate}")
    print(f"boltBearingStrengthYield: {boltBearingStrengthYield}")

    #Need to update axial force.
    #netAxialForce = 300
    netAxialForce = 16230
    #numberOfBolts = 12
    numberOfBolts = 20
    shearForcePerBolt = netAxialForce/numberOfBolts
    ultimateFOS = 2.0
    yieldFOS = 1.5

    if (boltShearStrength > boltBearingStrengthUltimate):
        print("Bearing Critical!")
        fittingFOS = 1.15
    else: 
        print("Sheer Critical!")
        fittingFOS = 2.0

    shearUltimateMOS = shearUMOS(boltShearStrength, shearForcePerBolt, ultimateFOS, fittingFOS)
    bearingUltimateMOS = bearingUMOS(boltBearingStrengthUltimate, shearForcePerBolt, ultimateFOS, fittingFOS)
    bearingYieldMOS = bearingYMOS(boltBearingStrengthYield, shearForcePerBolt, yieldFOS, fittingFOS)
    print(f"shearUltimateMOS: {shearUltimateMOS}")
    print(f"bearingUltimateMOS: {bearingUltimateMOS}")
    print(f"bearingYieldMOS: {bearingYieldMOS}")

main()