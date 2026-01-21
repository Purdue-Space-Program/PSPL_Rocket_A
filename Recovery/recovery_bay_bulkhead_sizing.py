pi = 3.14159
import sys
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *
import vehicle_parameters as v

def ShearStrength(material, diameter):
    return (material)*(pi*diameter*diameter*0.25)

def CalculateBearingStrength(material, diameter, thickness):
    bearing_strength = material * diameter * thickness
    return (bearing_strength)

def CalculateDesignLoad(load_per_bolt, FOS, fitting_factor):
    return (load_per_bolt*FOS*fitting_factor)

def CalculateMoS(maximum_allowable_load, limit_load, FOS, fitting_factor):
    MoS = (maximum_allowable_load/(limit_load*FOS*fitting_factor)) - 1
    return (MoS)

def main():

    fSU_AS_U = 180000.0
    F_bru_A_6061_T6 = 88_000 # [psi] from MMPDS-2019: https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/1934065665/Common+Material+Properties
    F_bry_A_6061_T6 = 58_000 # [psi] from MMPDS-2019: https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/1934065665/Common+Material+Properties
    plateThickness = 0.125
    ultimateFOS = 2.0
    yieldFOS = 1.5




    print("Recovery Bulkhead Bolted Joint")
    boltDiameter = 0.2075
    
    boltShearStrength = ShearStrength(fSU_AS_U, boltDiameter)
    boltBearingStrengthUltimate = CalculateBearingStrength(F_bru_A_6061_T6, boltDiameter, plateThickness)
    boltBearingStrengthYield = CalculateBearingStrength(F_bry_A_6061_T6, boltDiameter, plateThickness)
    print(f"\tboltShearStrength: {boltShearStrength}")
    print(f"\tboltBearingStrengthUltimate: {boltBearingStrengthUltimate:.2f}")
    print(f"\tboltBearingStrengthYield: {boltBearingStrengthYield:.2f}")
    
    netAxialForce = 3500
    numberOfBolts = 12
    shearForcePerBolt = netAxialForce/numberOfBolts

    if (boltShearStrength > boltBearingStrengthUltimate):
        print("\tBearing Critical!")
        fittingFOS = 1.15
    else: 
        print("\tSheer Critical!")
        fittingFOS = 2.0

    shearUltimateDesign = CalculateDesignLoad(shearForcePerBolt, ultimateFOS, fittingFOS)
    bearingUltimateDesign = CalculateDesignLoad(shearForcePerBolt, ultimateFOS, fittingFOS)
    bearingYieldDesign = CalculateDesignLoad(shearForcePerBolt, yieldFOS, fittingFOS)
    print(f"\tshearUltimateDesign: {shearUltimateDesign:.2f}")
    print(f"\tbearingUltimateDesign: {bearingUltimateDesign:.2f}")
    print(f"\tbearingYieldDesign: {bearingYieldDesign:.2f}")

    shearUltimateMOS = CalculateMoS(boltShearStrength, shearForcePerBolt, ultimateFOS, fittingFOS)
    bearingUltimateMOS = CalculateMoS(boltBearingStrengthUltimate, shearForcePerBolt, ultimateFOS, fittingFOS)
    bearingYieldMOS = CalculateMoS(boltBearingStrengthYield, shearForcePerBolt, yieldFOS, fittingFOS)
    print(f"\tshearUltimateMOS: {shearUltimateMOS:.2f}")
    print(f"\tbearingUltimateMOS: {bearingUltimateMOS:.2f}")
    print(f"\tbearingYieldMOS: {bearingYieldMOS:.2f}")

    
    
    
    
    print("Tank Bulkhead Bolted Joint")
    boltDiameter = 0.2614

    boltShearStrength = ShearStrength(fSU_AS_U, boltDiameter)
    boltBearingStrengthUltimate = CalculateBearingStrength(F_bru_A_6061_T6, boltDiameter, plateThickness)
    boltBearingStrengthYield = CalculateBearingStrength(F_bry_A_6061_T6, boltDiameter, plateThickness)
    print(f"\tboltShearStrength: {boltShearStrength}")
    print(f"\tboltBearingStrengthUltimate: {boltBearingStrengthUltimate:.2f}")
    print(f"\tboltBearingStrengthYield: {boltBearingStrengthYield:.2f}")
    
    netAxialForce = 16230
    numberOfBolts = 16
    shearForcePerBolt = netAxialForce/numberOfBolts

    if (boltShearStrength > boltBearingStrengthUltimate):
        print("\tBearing Critical!")
        fittingFOS = 1.15
    else: 
        print("\tSheer Critical!")
        fittingFOS = 2.0

    shearUltimateDesign = CalculateDesignLoad(shearForcePerBolt, ultimateFOS, fittingFOS)
    bearingUltimateDesign = CalculateDesignLoad(shearForcePerBolt, ultimateFOS, fittingFOS)
    bearingYieldDesign = CalculateDesignLoad(shearForcePerBolt, yieldFOS, fittingFOS)
    print(f"\tshearUltimateDesign: {shearUltimateDesign:.2f}")
    print(f"\tbearingUltimateDesign: {bearingUltimateDesign:.2f}")
    print(f"\tbearingYieldDesign: {bearingYieldDesign:.2f}")
    
    shearUltimateMOS = CalculateMoS(boltShearStrength, shearForcePerBolt, ultimateFOS, fittingFOS)
    bearingUltimateMOS = CalculateMoS(boltBearingStrengthUltimate, shearForcePerBolt, ultimateFOS, fittingFOS)
    bearingYieldMOS = CalculateMoS(boltBearingStrengthYield, shearForcePerBolt, yieldFOS, fittingFOS)
    print(f"\tshearUltimateMOS: {shearUltimateMOS:.2f}")
    print(f"\tbearingUltimateMOS: {bearingUltimateMOS:.2f}")
    print(f"\tbearingYieldMOS: {bearingYieldMOS:.2f}")

main()