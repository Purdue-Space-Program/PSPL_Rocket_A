import numpy as np
import sys
import os
from dataclasses import dataclass
import matplotlib.pyplot as plt
import copy
import scipy.io as sio
from pathlib import Path

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
import vehicle_parameters_functions
import vehicle_parameters
import vehicle_main
import print_filter


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

def CalculateCircleAreaWithDiameter(diameter):
    circle_area = np.pi*((diameter/2)**2)
    return (circle_area)

def CalculateMaximumAllowableBoltShearLoad(material_strength, 
                                           bolt_diameter):
    maximum_allowable_shear_load = material_strength * CalculateCircleAreaWithDiameter(bolt_diameter)
    return (maximum_allowable_shear_load)

def CalculateMaximumAllowableBearingLoad(joint_member_material_strength, 
                                         bolt_hole_diameter, 
                                         plate_thickness):    
    maximum_allowable_bearing_load = joint_member_material_strength * bolt_hole_diameter * plate_thickness
    return (maximum_allowable_bearing_load)

def CalculateMoS(maximum_allowable_load, 
                 limit_load, 
                 FoS, 
                 fitting_factor):
    MoS = (maximum_allowable_load/(limit_load*FoS*fitting_factor)) - 1
    return (MoS)


def Calculate_Shear_Bolted_Joint(bolt_thread_size, 
                                 bolt_material, 
                                 number_of_bolts, 
                                 shear_limit_load, 
                                 joint_member_1_material, 
                                 joint_member_1_thickness, 
                                 joint_member_1_E_d_ratio, 
                                 joint_member_1_shear_joint_type, 
                                 yield_FoS, 
                                 ultimate_FoS
                                 ):
    print(f"\tBolt Name: {bolt_thread_size} UNF")
    print(f"\tBolt Material: {bolt_material}")

    print(f"\tJoint Member 1 Material: {joint_member_1_material}")
    print(f"\tJoint Member 1 Thickness: {joint_member_1_thickness * c.M2IN} in")

    print(f"\tShear Joint Type: {joint_member_1_shear_joint_type}")
    match joint_member_1_shear_joint_type:
        case "Single":
            pass
        case "Double":
            shear_limit_load /= 2
        case _:
            raise ValueError(f"Invalid Shear Joint Type: {joint_member_1_shear_joint_type}")

    # material properties
    F_su_316_Stainless_Steel = 50_000 * c.PSI2PA # [psi] S-basis from MMPDS-2025
    F_su_Alloy_Steel = 153_000*0.5 * c.PSI2PA # [psi] out of my ass

    # from MMPDS-2019: https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/1934065665/Common+Material+Properties
    F_bru_Aluminum_6061_T6_E_d_one_point_five = 67_000 * c.PSI2PA # [psi]
    F_bru_Aluminum_6061_T6_E_d_two =            88_000 * c.PSI2PA # [psi]
    F_bry_Aluminum_6061_T6_E_d_one_point_five = 50_000 * c.PSI2PA # [psi]
    F_bry_Aluminum_6061_T6_E_d_two =            58_000 * c.PSI2PA # [psi]

    F_bru_Aluminum_7055_T74511_E_d_one_point_five = 115_000 * c.PSI2PA # [psi]
    F_bru_Aluminum_7055_T74511_E_d_two =            151_000 * c.PSI2PA # [psi]
    F_bry_Aluminum_7055_T74511_E_d_one_point_five = 96_000 * c.PSI2PA # [psi]
    F_bry_Aluminum_7055_T74511_E_d_two =            114_000 * c.PSI2PA # [psi]


    # F_bru_Aluminum_6063_T5_E_d_one_point_five = ???_000 * c.PSI2PA # [psi]
    F_bru_Aluminum_6063_T5_E_d_two = 46_000 * c.PSI2PA # [psi]
    # F_bry_Aluminum_6063_T5_E_d_one_point_five = ???_000 * c.PSI2PA # [psi]
    F_bry_Aluminum_6063_T5_E_d_two = 25_600 * c.PSI2PA # [psi]

    F_bru_304_Stainless_Steel_E_d_two = 262_000 * c.PSI2PA
    F_bry_304_Stainless_Steel_E_d_two = 123_000 * c.PSI2PA


    # source: https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/1838153742/magic+numbers
    match bolt_thread_size:
        case "#8":
            bolt_minor_diameter = 0.1299 * c.IN2M
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
                case "Aluminum 7055-T74511":
                    joint_member_F_bry = F_bry_Aluminum_7055_T74511_E_d_one_point_five
                    joint_member_F_bru = F_bru_Aluminum_7055_T74511_E_d_one_point_five
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
                case "Aluminum 7055-T74511":
                    joint_member_F_bry = F_bry_Aluminum_7055_T74511_E_d_two
                    joint_member_F_bru = F_bru_Aluminum_7055_T74511_E_d_two
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


    limit_shear_load_per_bolt = shear_limit_load/number_of_bolts
    print(f"limit_shear_load_per_bolt: {limit_shear_load_per_bolt * c.N2LBF:.2f} LBF")


    bolt_maximum_allowable_shear_ultimate_load = CalculateMaximumAllowableBoltShearLoad(bolt_F_su, bolt_minor_diameter)
    print(f"\tbolt_maximum_allowable_shear_ultimate_load: {bolt_maximum_allowable_shear_ultimate_load * c.N2LBF:.2f} LBF")

    # for member in joint_members
    #

    joint_member_1_maximum_allowable_bearing_yield_load = CalculateMaximumAllowableBearingLoad(joint_member_F_bry, bolt_hole_diameter, joint_member_1_thickness)
    joint_member_1_maximum_allowable_bearing_ultimate_load = CalculateMaximumAllowableBearingLoad(joint_member_F_bru, bolt_hole_diameter, joint_member_1_thickness)
    print(f"\tjoint_member_1_maximum_allowable_bearing_yield_load: {joint_member_1_maximum_allowable_bearing_yield_load * c.N2LBF:.2f} LBF, {joint_member_1_maximum_allowable_bearing_yield_load:.2f} N")
    print(f"\tjoint_member_1_maximum_allowable_bearing_ultimate_load: {joint_member_1_maximum_allowable_bearing_ultimate_load * c.N2LBF:.2f} LBF, {joint_member_1_maximum_allowable_bearing_ultimate_load:.2f} N")


    initial_fitting_factor = 1 # since we dont know if the joint is shear or bearing critical yet (idk if this is the right way to do it tbh)

    # bolt_shear_yield_MoS ????????????????????????? i have not found a shear yield strength value for any material so i guess it doesn't exist...
    bolt_shear_ultimate_MoS = CalculateMoS(maximum_allowable_load=bolt_maximum_allowable_shear_ultimate_load,
                                           limit_load=limit_shear_load_per_bolt,
                                           FoS=ultimate_FoS,
                                           fitting_factor=initial_fitting_factor)
    clamped_material_bearing_ultimate_MoS = CalculateMoS(maximum_allowable_load=joint_member_1_maximum_allowable_bearing_ultimate_load,
                                                         limit_load=limit_shear_load_per_bolt,
                                                         FoS=ultimate_FoS,
                                                         fitting_factor=initial_fitting_factor)

    if (bolt_maximum_allowable_shear_ultimate_load > joint_member_1_maximum_allowable_bearing_ultimate_load):
        fitting_factor = 1.15
        print(f"\tBearing Critical! 😄, Fitting Factor: {fitting_factor}")
    elif (bolt_maximum_allowable_shear_ultimate_load <= joint_member_1_maximum_allowable_bearing_ultimate_load):
        fitting_factor = 2.0
        print(f"\tShear Critical! 😢, Fitting Factor: {fitting_factor}")
    else:
        raise ValueError("what")


    bolt_shear_ultimate_MoS = CalculateMoS(bolt_maximum_allowable_shear_ultimate_load, limit_shear_load_per_bolt, ultimate_FoS, fitting_factor)
    clamped_material_bearing_ultimate_MoS = CalculateMoS(joint_member_1_maximum_allowable_bearing_ultimate_load, limit_shear_load_per_bolt, ultimate_FoS, fitting_factor)
    clamped_material_bearing_yield_MoS = CalculateMoS(joint_member_1_maximum_allowable_bearing_yield_load, limit_shear_load_per_bolt, yield_FoS, fitting_factor)
    # print(f"\tdesign_shear_load_per_bolt ULTIMATE: {limit_shear_load_per_bolt*(parameters.ultimate_FoS*fitting_factor) * c.N2LBF}")

    print(f"\t{MoS_Text_Color(bolt_shear_ultimate_MoS)}Bolt shear ultimate MoS: {bolt_shear_ultimate_MoS:.3f}{ENDC}", i_am_a_margin = True)
    print(f"\t{MoS_Text_Color(clamped_material_bearing_ultimate_MoS)}Clamped material bearing ultimate MoS: {clamped_material_bearing_ultimate_MoS:.3f}{ENDC}", i_am_a_margin = True)
    print(f"\t{MoS_Text_Color(clamped_material_bearing_yield_MoS)}Clamped material bearing yield MoS: {clamped_material_bearing_yield_MoS:.3f}{ENDC}", i_am_a_margin = True)

    print("", i_am_a_title = True)

    return (bolt_shear_ultimate_MoS, clamped_material_bearing_ultimate_MoS, clamped_material_bearing_yield_MoS)

@dataclass
class JointMember:
    material: str
    thickness: float
    E_d_ratio: float
    shear_joint_type: str

@dataclass
class ShearBoltedJoint:
    bolt_material: str
    bolt_thread_size: str
    number_of_bolts: float
    shear_limit_load: float
    joint_member_1: JointMember
    yield_FoS: float
    ultimate_FoS: float

    def Calculate_Shear_Bolted_Joint(self):
        return Calculate_Shear_Bolted_Joint(bolt_material =  self.bolt_material,
                                    bolt_thread_size = self.bolt_thread_size,
                                    number_of_bolts = self.number_of_bolts,
                                    shear_limit_load = self.shear_limit_load,
                                    joint_member_1_material = self.joint_member_1.material ,
                                    joint_member_1_thickness = self.joint_member_1.thickness,
                                    joint_member_1_E_d_ratio = self.joint_member_1.E_d_ratio,
                                    joint_member_1_shear_joint_type = self.joint_member_1.shear_joint_type,
                                    yield_FoS = self.yield_FoS,
                                    ultimate_FoS= self.ultimate_FoS,
                                    )

def Calculate_Shear_Bolted_Joints(parameters):

    upper_strut = JointMember(material = "Aluminum 6063-T52",
                            thickness = 0.125 * c.IN2M,
                            E_d_ratio = 2.0,
                            shear_joint_type = "Double"
    )
    mid_strut = copy.deepcopy(upper_strut)
    lower_strut = copy.deepcopy(upper_strut)
    lower_strut.shear_joint_type = "Single"

    tank_wall = JointMember(material = "Aluminum 6061-T6",
                            thickness = 0.125 * c.IN2M,
                            E_d_ratio = 1.5,
                            shear_joint_type = "Single"
                           )

    recovery_bay = copy.deepcopy(tank_wall)


    print("-------------Tank Wall to Bulkhead Bolted Joint-------------", i_am_a_title=True)
    bulkhead_area = CalculateCircleAreaWithDiameter(parameters.tank_inner_diameter)
    bulkhead_nominal_pressure_blowoff_limit_load = (parameters.nominal_tank_pressure * bulkhead_area) * parameters.proof_factor
    bulkhead_largest_possible_pressure_blowoff_limit_load = (parameters.largest_possible_tank_pressure * bulkhead_area) * parameters.proof_factor
    print(f"\tBulkhead blowoff load: {bulkhead_largest_possible_pressure_blowoff_limit_load:.2f} N, {bulkhead_largest_possible_pressure_blowoff_limit_load * c.N2LBF :.2f} LBF")
    # print(f"\tOxygen tank max load: {parameters.oxygen_tank_max_load:.2f} N, {parameters.oxygen_tank_max_load * c.N2LBF :.2f} LBF")

    bulkhead_max_limit_load = bulkhead_largest_possible_pressure_blowoff_limit_load
    print("\tbulkhead_max_load: bulkhead_blowoff_load")

    # if bulkhead_blowoff_limit_load > parameters.oxygen_tank_max_load:
    #     bulkhead_max_limit_load = bulkhead_blowoff_limit_load
    #     print("\tbulkhead_max_load: bulkhead_blowoff_load")
    # else:
    #     bulkhead_max_limit_load = parameters.oxygen_tank_max_load
    #     print("\tbulkhead_max_load: parameters.oxygen_tank_max_load")

    tank_wall_to_bulkhead_joint = ShearBoltedJoint(bolt_material = "Alloy Steel",
                                                   bolt_thread_size = "5/16\"",
                                                   number_of_bolts = 30,
                                                   shear_limit_load = bulkhead_max_limit_load,
                                                   joint_member_1 = tank_wall,
                                                   yield_FoS = parameters.yield_FoS,
                                                   ultimate_FoS = parameters.ultimate_FoS,
                                                  )
    tank_wall_to_bulkhead_joint.Calculate_Shear_Bolted_Joint()


    print("------------- Tank Bulkhead to Upper Strut Bolted Joint -------------", i_am_a_title=True)
    tank_bulkhead_to_upper_strut_joint = ShearBoltedJoint(bolt_material = "316 Stainless Steel",
                                                          bolt_thread_size = "#10",
                                                          number_of_bolts = 2,
                                                          shear_limit_load = parameters.upper_strut_max_load,
                                                          joint_member_1 = upper_strut,
                                                          yield_FoS = parameters.yield_FoS,
                                                          ultimate_FoS = parameters.ultimate_FoS,
                                                         )
    tank_bulkhead_to_upper_strut_joint.Calculate_Shear_Bolted_Joint()

    print("------------- Tank Bulkhead to Mid Strut Bolted Joint -------------", i_am_a_title=True)
    tank_bulkhead_to_mid_strut_joint = copy.deepcopy(tank_bulkhead_to_upper_strut_joint)
    tank_bulkhead_to_mid_strut_joint.joint_member_1 = mid_strut
    tank_bulkhead_to_mid_strut_joint.shear_limit_load = parameters.mid_strut_max_load
    tank_bulkhead_to_mid_strut_joint.Calculate_Shear_Bolted_Joint()

    print("------------- Tank Bulkhead to Lower Strut Bolted Joint -------------", i_am_a_title=True)
    tank_bulkhead_to_lower_strut_joint = copy.deepcopy(tank_bulkhead_to_upper_strut_joint)
    tank_bulkhead_to_lower_strut_joint.shear_limit_load = parameters.lower_strut_max_load
    tank_bulkhead_to_lower_strut_joint.joint_member_1 = lower_strut
    tank_bulkhead_to_lower_strut_joint.Calculate_Shear_Bolted_Joint()


    print("------------- Recovery Bulkhead Bolted Joint -------------", i_am_a_title=True)
    injector_upper_half_to_fin_can_strut = ShearBoltedJoint(bolt_material = "316 Stainless Steel",
                                                            bolt_thread_size = "1/4",
                                                            number_of_bolts = 12,
                                                            shear_limit_load = parameters.copv_tube_max_load,
                                                            joint_member_1 = recovery_bay,
                                                            yield_FoS = parameters.yield_FoS,
                                                            ultimate_FoS = parameters.ultimate_FoS,
                                                            )
    injector_upper_half_to_fin_can_strut.Calculate_Shear_Bolted_Joint()



    # print("-------------Recovery Bay Connector Bolted Joint-------------\n", i_am_a_title=True)


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

    return(parameters)


def main(parameters):
    parameters = Calculate_Shear_Bolted_Joints(parameters)
    return(parameters)

if __name__ == "__main__":
    rerun_everything = False

    if rerun_everything:
        vehicle_main.vehicle_analysis()
    else:
        pass
        # an idea i have is that the alternative to running everything is using the last saved parameters and giving a warning say it might not be as accurate/up to date
        # parameters = csv_to_dataclass(parameters_csv_filepath)

        parameters, wet_mass_distribution, dry_mass_distribution = vehicle_parameters.main()
        with print_filter.context_manager(print_everything=True):
            main(parameters)    
