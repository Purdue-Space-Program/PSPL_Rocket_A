import numpy as np
import sys
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
from vehicle_parameters import parameters as parameters


def CalculateCircleAreaWithDiameter(diameter):
    circle_area = np.pi*((diameter/2)**2)
    return (circle_area)

def CalculateMaximumAllowableBoltShearLoad(material_strength, bolt_diameter):
    maximum_allowable_shear_load = material_strength * CalculateCircleAreaWithDiameter(bolt_diameter)
    return (maximum_allowable_shear_load)

def CalculateMaximumAllowableBearingLoad(material_strength, bolt_hole_diameter, plate_thickness):
    maximum_allowable_bearing_load = material_strength * bolt_hole_diameter * plate_thickness
    return (maximum_allowable_bearing_load)

def CalculateMoS(maximum_allowable_load, limit_load, FOS, fitting_factor):
    MoS = (maximum_allowable_load/(limit_load*FOS*fitting_factor)) - 1
    return (MoS)

def main():
    print("Recovery Bay Connector Bolted Joint")
    
    # bolt dimensions
    bolt_name = "5/16"
    number_of_bolts = 18
    E_d_ratio = 1.5
    
    ultimate_FOS = 2.0
    yield_FOS = 1.5    
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
    
    
    # source: https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/1838153742/magic+numbers
    match bolt_name:
        case "#10":
            bolt_minor_diameter = 0.1517 * c.IN2M
            bolt_hole_clearance_diameter_tight = 0.1960 * c.IN2M
            bolt_hole_clearance_diameter_loose = 0.2010 * c.IN2M
        case "1/4":
            bolt_minor_diameter = 0.2075 * c.IN2M
            bolt_hole_clearance_diameter_tight = 0.2570 * c.IN2M
        case "5/16":
            bolt_minor_diameter = 0.2614 * c.IN2M
        case "3/8":
            bolt_minor_diameter = 0.32 * c.IN2M
        case _:
            raise ValueError("balls")
    
    bolt_hole_diameter = bolt_minor_diameter

    
    print(f"Bolt Name: {bolt_name} UNF")
    bolt_maximum_allowable_shear_ultimate_load = CalculateMaximumAllowableBoltShearLoad(F_su_SS_316, bolt_minor_diameter)
    print(f"\tboltShearStrength: {bolt_maximum_allowable_shear_ultimate_load:.2f} ")

    
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
    net_axial_force = blowoff_load
    shear_force_per_bolt = net_axial_force/number_of_bolts

        
    initial_fitting_factor = 1.15 # since we dont know if the joint is shear or bearing critical yet (idk if this is the right way to do it tbh)
    
    # bolt_shear_yield_MoS ?????????????????????????
    bolt_shear_ultimate_MoS = CalculateMoS(bolt_maximum_allowable_shear_ultimate_load, shear_force_per_bolt, ultimate_FOS, initial_fitting_factor)
    tank_wall_bearing_ultimate_MoS = CalculateMoS(tank_wall_maximum_allowable_bearing_ultimate_load, shear_force_per_bolt, ultimate_FOS, initial_fitting_factor)
    
    if (bolt_maximum_allowable_shear_ultimate_load > tank_wall_maximum_allowable_bearing_ultimate_load):
        print("\tBearing Critical! 😄")
        fitting_factor = 1.15
    elif (bolt_maximum_allowable_shear_ultimate_load <= tank_wall_maximum_allowable_bearing_ultimate_load):
        print("\tShear Critical!")
        fitting_factor = 2.0
    else:
        raise ValueError("what")
        

    bolt_shear_ultimate_MoS = CalculateMoS(bolt_maximum_allowable_shear_ultimate_load, shear_force_per_bolt, ultimate_FOS, fitting_factor)
    tank_wall_bearing_ultimate_MoS = CalculateMoS(tank_wall_maximum_allowable_bearing_ultimate_load, shear_force_per_bolt, ultimate_FOS, fitting_factor)
    bolt_shear_ultimate_MoS = CalculateMoS(bolt_maximum_allowable_shear_ultimate_load, shear_force_per_bolt, ultimate_FOS, fitting_factor)
    tank_wall_bearing_ultimate_MoS = CalculateMoS(tank_wall_maximum_allowable_bearing_ultimate_load, shear_force_per_bolt, ultimate_FOS, fitting_factor)
    bearingYieldMOS = CalculateMoS(tank_wall_maximum_allowable_bearing_yield_load, shear_force_per_bolt, yield_FOS, fitting_factor)
    
    print(f"\tshearUltimateMOS: {bolt_shear_ultimate_MoS:.3f}")
    print(f"\tbearingUltimateMOS: {tank_wall_bearing_ultimate_MoS:.3f}")
    print(f"\tbearingYieldMOS: {bearingYieldMOS:.3f}")

main()