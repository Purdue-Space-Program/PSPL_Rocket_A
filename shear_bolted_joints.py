import numpy as np
import sys
import os
from dataclasses import dataclass, fields, field
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
from vehicle_parameters import parameters
import print_filter

def CalculateCircleAreaWithDiameter(diameter):
    circle_area = np.pi*((diameter/2)**2)
    return (circle_area)

def CalculateMaximumAllowableBoltShearLoad(material_strength, bolt_diameter):
    maximum_allowable_shear_load = material_strength * CalculateCircleAreaWithDiameter(bolt_diameter)
    return (maximum_allowable_shear_load)

def CalculateMaximumAllowableBearingLoad(joint_member_material_strength, bolt_hole_diameter, plate_thickness):
    maximum_allowable_bearing_load = joint_member_material_strength * bolt_hole_diameter * plate_thickness
    return (maximum_allowable_bearing_load)

def CalculateMoS(maximum_allowable_load, limit_load, FOS, fitting_factor):
    MoS = (maximum_allowable_load/(limit_load*FOS*fitting_factor)) - 1
    return (MoS)


def Calculate_Shear_Bolted_Joint(bolt_thread_size, bolt_material, number_of_bolts, joint_member_1_material, joint_member_1_thickness, joint_member_1_shear_limit_load, joint_member_1_E_d_ratio, joint_member_1_shear_joint_type):
    print(f"\tBolt Name: {bolt_thread_size} UNF")
    print(f"\tBolt Material: {bolt_material}")
    
    print(f"\tJoint Member 1 Material: {joint_member_1_material}")
    print(f"\tJoint Member 1 Thickness: {joint_member_1_thickness}")
    
    print(f"\tShear Joint Type: {joint_member_1_shear_joint_type}")
    match joint_member_1_shear_joint_type:
        case "Single":
            pass
        case "Double":
            joint_member_1_shear_limit_load /= 2
        case _:
            raise ValueError(f"Invalid Shear Joint Type: {joint_member_1_shear_joint_type}")

    # material properties
    F_su_316_Stainless_Steel = 50_000 * c.PSI2PA # [psi] S-basis from MMPDS-2019
    F_su_Alloy_Steel = 153_000*0.6 * c.PSI2PA # [psi] S-basis from MMPDS-2019

    # from MMPDS-2019: https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/1934065665/Common+Material+Properties
    F_bry_Aluminum_6061_T6_E_d_one_point_five = 50_000 * c.PSI2PA # [psi] 
    F_bru_Aluminum_6061_T6_E_d_one_point_five = 67_000 * c.PSI2PA # [psi]
    F_bry_Aluminum_6061_T6_E_d_two = 58_000 * c.PSI2PA # [psi] 
    F_bru_Aluminum_6061_T6_E_d_two = 88_000 * c.PSI2PA # [psi] 
    
    # F_bry_Aluminum_6063_T5_E_d_one_point_five = ???_000 * c.PSI2PA # [psi] 
    # F_bru_Aluminum_6063_T5_E_d_one_point_five = ???_000 * c.PSI2PA # [psi]
    F_bry_Aluminum_6063_T5_E_d_two = 25_600 * c.PSI2PA # [psi] 
    F_bru_Aluminum_6063_T5_E_d_two = 46_000 * c.PSI2PA # [psi] 

    F_bry_304_Stainless_Steel_E_d_two = 123_000 * c.PSI2PA
    F_bru_304_Stainless_Steel_E_d_two = 262_000 * c.PSI2PA

    
    # source: https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/1838153742/magic+numbers
    match bolt_thread_size:
        case "#10":
            bolt_minor_diameter = 0.1517 * c.IN2M
            bolt_hole_clearance_diameter_tight = 0.1960 * c.IN2M
            bolt_hole_clearance_diameter_loose = 0.2010 * c.IN2M
        case "1/4":
            bolt_minor_diameter = 0.2075 * c.IN2M
            bolt_hole_clearance_diameter_tight = 0.2570 * c.IN2M
        case "5/16\"":
            bolt_minor_diameter = 0.2614 * c.IN2M
        case "3/8":
            bolt_minor_diameter = 0.32 * c.IN2M
        case _:
            raise ValueError("balls")
    
    bolt_hole_diameter = bolt_minor_diameter
    
    
    
    match joint_member_1_E_d_ratio:
        case 1.5:
            match joint_member_1_material:
                case "Aluminum 6061-T6":
                    joint_member_F_bry = F_bry_Aluminum_6061_T6_E_d_one_point_five
                    joint_member_F_bru = F_bru_Aluminum_6061_T6_E_d_one_point_five
                # case "Aluminum 6063-T52":
                    # joint_member_F_bry = F_bry_Aluminum_6063_T5_E_d_one_point_five
                    # joint_member_F_bru = F_bru_Aluminum_6063_T5_E_d_one_point_five
                case _:
                    raise ValueError("cock and ball torque")
        case 2.0:
            match joint_member_1_material:
                case "Aluminum 6061-T6":
                    joint_member_F_bry = F_bry_Aluminum_6061_T6_E_d_two
                    joint_member_F_bru = F_bru_Aluminum_6061_T6_E_d_two
                case "Aluminum 6063-T52":
                    joint_member_F_bry = F_bry_Aluminum_6063_T5_E_d_two
                    joint_member_F_bru = F_bru_Aluminum_6063_T5_E_d_two
                case "304 Stainless Steel":
                    joint_member_F_bry = F_bry_304_Stainless_Steel_E_d_two
                    joint_member_F_bru = F_bru_304_Stainless_Steel_E_d_two
                case _:
                    raise ValueError("cock and ball torque")
        case _:
            raise ValueError("cock and ball torque")
    
    match bolt_material:
        case "316 Stainless Steel":
            bolt_F_su = F_su_316_Stainless_Steel
        case "Alloy Steel":
            bolt_F_su = F_su_Alloy_Steel
        case _:
            raise ValueError("cock and ball torque")
            
    
    bolt_maximum_allowable_shear_ultimate_load = CalculateMaximumAllowableBoltShearLoad(bolt_F_su, bolt_minor_diameter)
    print(f"\tBolt shear strength: {bolt_maximum_allowable_shear_ultimate_load:.2f} ")
    
    
    joint_member_1_maximum_allowable_bearing_yield_load = CalculateMaximumAllowableBearingLoad(joint_member_F_bry, bolt_hole_diameter, joint_member_1_thickness)
    joint_member_1_maximum_allowable_bearing_ultimate_load = CalculateMaximumAllowableBearingLoad(joint_member_F_bru, bolt_hole_diameter, joint_member_1_thickness)
    print(f"\tjoint_member_1_maximum_allowable_bearing_yield_load: {joint_member_1_maximum_allowable_bearing_yield_load * c.N2LBF:.2f} LBF, {joint_member_1_maximum_allowable_bearing_yield_load:.2f} N")
    print(f"\tjoint_member_1_maximum_allowable_bearing_ultimate_load: {joint_member_1_maximum_allowable_bearing_ultimate_load * c.N2LBF:.2f} LBF, {joint_member_1_maximum_allowable_bearing_ultimate_load:.2f} N")

    limit_shear_load_per_bolt = joint_member_1_shear_limit_load/number_of_bolts
    
    initial_fitting_factor = 1 # since we dont know if the joint is shear or bearing critical yet (idk if this is the right way to do it tbh)
    
    # bolt_shear_yield_MoS ????????????????????????? i have not found a shear yield strength value for any material so i guess it doesnt exist...
    bolt_shear_ultimate_MoS = CalculateMoS(bolt_maximum_allowable_shear_ultimate_load, limit_shear_load_per_bolt, parameters.ultimate_FoS, initial_fitting_factor)
    clamped_material_bearing_ultimate_MoS = CalculateMoS(joint_member_1_maximum_allowable_bearing_ultimate_load, limit_shear_load_per_bolt, parameters.ultimate_FoS, initial_fitting_factor)
    
    if (bolt_maximum_allowable_shear_ultimate_load > joint_member_1_maximum_allowable_bearing_ultimate_load):
        fitting_factor = 1.15
        print(f"\tBearing Critical! 😄, Fitting Factor: {fitting_factor}")
    elif (bolt_maximum_allowable_shear_ultimate_load <= joint_member_1_maximum_allowable_bearing_ultimate_load):
        fitting_factor = 2.0
        print(f"\tShear Critical! 😢, Fitting Factor: {fitting_factor}")
    else:
        raise ValueError("what")
    
    
    bolt_shear_ultimate_MoS = CalculateMoS(bolt_maximum_allowable_shear_ultimate_load, limit_shear_load_per_bolt, parameters.ultimate_FoS, fitting_factor)
    clamped_material_bearing_yield_MoS = CalculateMoS(joint_member_1_maximum_allowable_bearing_yield_load, limit_shear_load_per_bolt, parameters.yield_FoS, fitting_factor)
    clamped_material_bearing_ultimate_MoS = CalculateMoS(joint_member_1_maximum_allowable_bearing_ultimate_load, limit_shear_load_per_bolt, parameters.ultimate_FoS, fitting_factor)
    
    
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    ENDC = '\033[0m'
    GREEN = '\033[92m'
    
    def MoS_Text_Color(MoS):
        if MoS >= 0:
            text_color = GREEN
        else:
            text_color = RED
        return(text_color)
    
    print(f"\t{MoS_Text_Color(bolt_shear_ultimate_MoS)}Bolt shear ultimate MoS: {bolt_shear_ultimate_MoS:.3f}{ENDC}", i_am_a_margin = True)
    print(f"\t{MoS_Text_Color(clamped_material_bearing_ultimate_MoS)}Clamped material bearing ultimate MoS: {clamped_material_bearing_ultimate_MoS:.3f}{ENDC}", i_am_a_margin = True)
    print(f"\t{MoS_Text_Color(clamped_material_bearing_yield_MoS)}Clamped material bearing yield MoS: {clamped_material_bearing_yield_MoS:.3f}{ENDC}", i_am_a_margin = True)

    print("", i_am_a_title = True)
    
    return (bolt_shear_ultimate_MoS, clamped_material_bearing_ultimate_MoS, clamped_material_bearing_yield_MoS)



@dataclass
class JointMember:
    # General Parameters
    material: str
    thickness: float
    shear_limit_load: float
    E_d_ratio: float
    shear_joint_type: str


@dataclass
class ShearBoltedJoint:
    # General Parameters
    bolt_material: str
    bolt_thread_size: str
    number_of_bolts: float
    joint_member_1: JointMember
    # joint_member_1_material: str
    # joint_member_1_thickness: float
    # joint_member_1_shear_limit_load: float
    # E_d_ratio: float
    # shear_joint_type: str
        
    def Calculate_Shear_Bolted_Joint(self):
        return Calculate_Shear_Bolted_Joint(bolt_material =  self.bolt_material, 
                                    bolt_thread_size = self.bolt_thread_size, 
                                    number_of_bolts = self.number_of_bolts,
                                    joint_member_1_material = self.joint_member_1.material ,
                                    joint_member_1_thickness = self.joint_member_1.thickness,
                                    joint_member_1_shear_limit_load = self.joint_member_1.shear_limit_load,        
                                    joint_member_1_E_d_ratio = self.joint_member_1.E_d_ratio,
                                    joint_member_1_shear_joint_type = self.joint_member_1.shear_joint_type)

def Calculate_Shear_Bolted_Joints():
    
    lower_strut = JointMember(material = "Aluminum 6063-T52",
                            thickness = 0.125 * c.IN2M,
                            E_d_ratio = 1.5,
                            shear_limit_load = 1 * c.LBF2N,
                            shear_joint_type = "Double"
    )
    mid_strut = JointMember(material = "Aluminum 6063-T52",
                            thickness = 0.125 * c.IN2M,
                            E_d_ratio = 2,
                            shear_limit_load = 991 * c.LBF2N,
                            shear_joint_type = "Double"
    )
    upper_strut = JointMember(material = "Aluminum 6063-T52",
                            thickness = 0.125 * c.IN2M,
                            E_d_ratio = 1.5,
                            shear_limit_load = 991 * c.LBF2N,
                            shear_joint_type = "Double"
    )
    
    tank_wall = JointMember(material = "Aluminum 6061-T6",
                            thickness = 0.125 * c.IN2M,
                            E_d_ratio = 1.5,
                            shear_limit_load = 991 * c.LBF2N,
                            shear_joint_type = "Single"
    )
    
    
    # print("-------------Tank Wall to Bulkhead Bolted Joint-------------", i_am_a_title=True)
    # bulkhead_area = CalculateCircleAreaWithDiameter(parameters.tank_inner_diameter)
    # bulkhead_blowoff_load = (parameters.tank_pressure * bulkhead_area) * parameters.proof_factor
    # print(f"\tBulkhead blowoff load: {bulkhead_blowoff_load:.2f} N, {bulkhead_blowoff_load * c.N2LBF :.2f} LBF")
        
    # tank_wall_to_bulkhead_joint = ShearBoltedJoint(bolt_material = "Alloy Steel", 
    #                                                     bolt_thread_size = "5/16\"", 
    #                                                     number_of_bolts = 18,
    #                                                     joint_member_1_material = "Aluminum 6061-T6",
    #                                                     joint_member_1_thickness = 0.125 * c.IN2M,
    #                                                     E_d_ratio = 1.5,
    #                                                     joint_member_1_shear_limit_load = bulkhead_blowoff_load,
    #                                                     shear_joint_type = "Single")
    # tank_wall_to_bulkhead_joint.Calculate_Shear_Bolted_Joint()


    
    print("------------- Tank Bulkhead to Mid Strut Bolted Joint -------------", i_am_a_title=True)
    tank_bulkhead_to_strut_joint = ShearBoltedJoint(bolt_material = "316 Stainless Steel", 
                                                        bolt_thread_size = "#10", 
                                                        number_of_bolts = 2,
                                                        joint_member_1 = mid_strut)
    tank_bulkhead_to_strut_joint.Calculate_Shear_Bolted_Joint()


    # print("-------------Launch Lug Bolted Joint-------------", i_am_a_title=True)
    # tank_bulkhead_to_strut_joint = ShearBoltedJoint(bolt_material = "316 Stainless Steel", 
    #                                                     bolt_thread_size = "#10", 
    #                                                     number_of_bolts = 2,
    #                                                     # joint_member_1 = launch_lug)
    #                                                     joint_member_1_material = "Aluminum 6061-T6", # https://www.speedymetals.com/pc-4676-8379-34-sq-wall-sq-tube-6063-t52-aluminum.aspx
    #                                                     joint_member_1_thickness = 0.5 * c.IN2M,
    #                                                     E_d_ratio = 2,
    #                                                     joint_member_1_shear_limit_load = 668 * c.LBF2N,
    #                                                     shear_joint_type = "Single")
    # tank_bulkhead_to_strut_joint.Calculate_Shear_Bolted_Joint()



    print("------------- Recovery Bulkhead Bolted Joint -------------", i_am_a_title=True)    
    injector_upper_half_to_fin_can_strut = ShearBoltedJoint(bolt_material = "316 Stainless Steel", 
                                                            bolt_thread_size = "1/4", 
                                                            number_of_bolts = 12,
                                                            joint_member_1 = tank_wall,
                                                            )
    injector_upper_half_to_fin_can_strut.Calculate_Shear_Bolted_Joint()
    
    
    
    # print("-------------Recovery Bay Connector Bolted Joint-------------\n", i_am_a_title=True)
    # # tank wall
    
    # tank_wall_to_bulkhead_bolted_joint_bolt_material = "316 Stainless Steel"
    # tank_wall_to_bulkhead_bolted_joint_joint_member_material = "Aluminum 6061-T6"

    # tank_wall_to_bulkhead_bolted_joint_bolt_name = "5/16"
    # tank_wall_to_bulkhead_bolted_joint_number_of_bolts = 18
    # tank_wall_to_bulkhead_bolted_joint_E_d_ratio = 1.5

    # tank_wall_thickness = 0.125 * c.IN2M

    # # bolt
    # tank_wall_to_bulkhead_bolted_joint_number_of_bolts = 18
    # tank_wall_to_bulkhead_bolted_joint_E_d_ratio = 1.5

    # CalculateShearJoint(tank_wall_to_bulkhead_bolted_joint_bolt_name, 
    #                     tank_wall_to_bulkhead_bolted_joint_bolt_material, 
    #                     tank_wall_to_bulkhead_bolted_joint_number_of_bolts, 
    #                     tank_wall_to_bulkhead_bolted_joint_joint_member_material, 
    #                     tank_wall_thickness, 
    #                     tank_wall_to_bulkhead_bolted_joint_E_d_ratio)
    


def main():
    with print_filter.context_manager(print_everything=False):
        Calculate_Shear_Bolted_Joints()


if __name__ == "__main__":
    main()
    