import numpy as np
import sys
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
from vehicle_parameters import parameters as parameters


def CalculateCircleAreaWithDiameter(diameter):
    circle_area = np.pi*(diameter/2)**2
    return (circle_area)

def CalculateMaximumAllowableBoltShearLoad(material_strength, bolt_diameter):
    maximum_allowable_shear_load = material_strength * CalculateCircleAreaWithDiameter(bolt_diameter)
    return (maximum_allowable_shear_load)

def CalculateMaximumAllowableBearingLoad(material_strength, bolt_hole_diameter, plate_thickness):
    maximum_allowable_bearing_load = material_strength * bolt_hole_diameter * plate_thickness
    return (maximum_allowable_bearing_load)

def CalculateMOS(maximum_allowable_load, limit_load, FOS, fitting_factor):
    MoS = (maximum_allowable_load/(limit_load*FOS*fitting_factor)) - 1
    return (MoS)

def main():
    print("Recovery Bay Connector Bolted Joint")
    
    # bolt dimensions
    bolt_name = "5/16"
    number_of_bolts = 18
    E_d_ratio = 1.5
    
    ultimateFOS = 2.0
    yieldFOS = 1.5    
    proof_factor = 1.5
    
    # material properties
    F_su_SS_316 = 66_000 * c.PSI2PA # [psi]
    F_su_ALLOY_STEEL = 180_000 * c.PSI2PA # [psi]

    
    
    match E_d_ratio:
        case 1.5:
            F_bry_A_6061_T6 = 50_000 * c.PSI2PA # [psi] from MMPDS-2019: https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/1934065665/Common+Material+Properties
            F_bru_A_6061_T6 = 67_000 * c.PSI2PA # [psi] from MMPDS-2019: https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/1934065665/Common+Material+Properties
        case 2.0:
            F_bry_A_6061_T6 = 58_000 * c.PSI2PA # [psi] from MMPDS-2019: https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/1934065665/Common+Material+Properties
            F_bru_A_6061_T6 = 88_000 * c.PSI2PA # [psi] from MMPDS-2019: https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/1934065665/Common+Material+Properties
        case _:
            raise ValueError("cock and ball torque")
    
    
    match bolt_name:
        case "1/4":
            bolt_minor_diameter = 0.2075 * c.IN2M # [in]
        case "5/16":
            bolt_minor_diameter = 0.2614 * c.IN2M # [in]
        case "3/8":
            bolt_minor_diameter = 0.32 * c.IN2M # [in]
        case _:
            raise ValueError("balls")
    
    bolt_hole_diameter = bolt_minor_diameter

    
    print(f"Bolt Name: {bolt_name} UNF")
    maximum_allowable_bolt_shear_ultimate_load = CalculateMaximumAllowableBoltShearLoad(F_su_SS_316, bolt_minor_diameter)
    print(f"\tboltShearStrength: {maximum_allowable_bolt_shear_ultimate_load:.2f} ")

    
    tank_wall_thickness = 0.125 * c.IN2M
    tank_wall_maximum_allowable_bearing_yield_load = CalculateMaximumAllowableBearingLoad(F_bry_A_6061_T6, bolt_hole_diameter, tank_wall_thickness)
    tank_wall_maximum_allowable_bearing_ultimate_load = CalculateMaximumAllowableBearingLoad(F_bru_A_6061_T6, bolt_hole_diameter, tank_wall_thickness)
    print(f"\tank_wall_maximum_allowable_bearing_yield_load: {tank_wall_maximum_allowable_bearing_yield_load:.2f}")
    print(f"\tank_wall_maximum_allowable_bearing_ultimate_load: {tank_wall_maximum_allowable_bearing_ultimate_load * c.N2LBF:.2f} LBF, {tank_wall_maximum_allowable_bearing_ultimate_load:.2f} N")

    bulkhead_area = CalculateCircleAreaWithDiameter(parameters.tank_inner_diameter)
    blowoff_load = parameters.tank_pressure * bulkhead_area * proof_factor
    print(f"\tblowoff_load: {blowoff_load:.2f} N, {blowoff_load * c.N2LBF :.2f} LBF")
    
    #Need to update axial force.
    #netAxialForce = 300
    netAxialForce = blowoff_load
    #numberOfBolts = 12
    shearForcePerBolt = netAxialForce/number_of_bolts

        
    initial_fitting_factor = 1.15 # since we dont know if the joint is shear or bearing critical yet
    
    shearUltimateMOS = CalculateMOS(maximum_allowable_bolt_shear_ultimate_load, shearForcePerBolt, ultimateFOS, initial_fitting_factor)
    bearingUltimateMOS = CalculateMOS(tank_wall_maximum_allowable_bearing_ultimate_load, shearForcePerBolt, ultimateFOS, initial_fitting_factor)
    
    if (maximum_allowable_bolt_shear_ultimate_load > tank_wall_maximum_allowable_bearing_ultimate_load):
        print("\tBearing Critical! 😄")
        fitting_factor = 1.15
    elif (maximum_allowable_bolt_shear_ultimate_load <= tank_wall_maximum_allowable_bearing_ultimate_load):
        print("\tSheer Critical!")
        fitting_factor = 2.0
    else:
        raise ValueError("what")
        

    shearUltimateMOS = CalculateMOS(maximum_allowable_bolt_shear_ultimate_load, shearForcePerBolt, ultimateFOS, fitting_factor)
    bearingUltimateMOS = CalculateMOS(tank_wall_maximum_allowable_bearing_ultimate_load, shearForcePerBolt, ultimateFOS, fitting_factor)
    bearingYieldMOS = CalculateMOS(tank_wall_maximum_allowable_bearing_yield_load, shearForcePerBolt, yieldFOS, fitting_factor)
    print(f"\tshearUltimateMOS: {shearUltimateMOS:.3f}")
    print(f"\tbearingUltimateMOS: {bearingUltimateMOS:.3f}")
    print(f"\tbearingYieldMOS: {bearingYieldMOS:.3f}")

main()