import sys
import os
import matplotlib.pyplot as plt
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *
import vehicle_parameters


def tensile_from_area(diameter, pressure):
    return ((diameter*diameter/4.0*np.pi) * pressure)

def tensile_from_o_ring(diameter, required_compression_per_linear_inch):
    return (diameter*np.pi) * required_compression_per_linear_inch

#def recovery_strut_force(mass, acceleration)

def proof(safety, strength):
    return safety*strength

def upper_preload(diameter, proof, preload_percent):
    area = np.pi * (diameter*diameter/4)
    return ((area*proof)*preload_percent)

def lower_preload(diameter, proof, preload_percent, torque_variation):
    return (upper_preload(diameter, proof, preload_percent)/(1+torque_variation))*(1-torque_variation)

def bolts(force, preload):
    return force/preload

def sum_forces(array, safety):
    net = 0
    for i in array:
        net = net + i
    return net * safety   

def diameter(ratio, bolt, nut, ID, wall):
    radius = ID/2 + wall + nut/2 + bolt*ratio
    return 2*radius

def MOS(calculated, actual):
    return (actual - calculated)/calculated










if __name__ == "__main__":
    

    parameters, wet_mass_distribution, dry_mass_distribution = vehicle_parameters.main()    
    proof_stress = 140000
    safety_factor = 1.4*1.15
    preload_percent = 0.75
    torque_variation = 0.25
    #Pathfinder
    bolt_major_diameter_main = 0.25 #0.19
    bolt_minor_diameter_main = 0.2075 #0.1528
    bolt_major_diameter_pintle = 0.125
    bolt_minor_diameter_pintle = 0.0979
    chamber_diameter = 5
    chamber_wall_thickness = 0.25
    chamber_pressure = 200
    #CMS
    #bolt_diameter_plates_major = 0.25
    #bolt_diameter_plates_minor = 0.2075
    #chamber_diameter = 5.05
    #chamber_pressure = 200

    nut_diameter_main = bolt_major_diameter_main*2.5
    ED_ratio = 1.5
    chamber_wall_thickness = 0.25


    force_manifold_pressure = tensile_from_area(chamber_diameter,chamber_pressure/0.8)
    print(f"force_injector_pressure: {force_manifold_pressure:.2f}")

    force_pintle_pressure = tensile_from_area(.98,chamber_pressure/0.8)
    print(f"force_pintle_pressure: {force_pintle_pressure:.2f}")


    force_outer_o_ring = tensile_from_o_ring(5.44301,40)
    force_film_o_ring = tensile_from_o_ring(4.193,40)
    force_PT_o_ring = tensile_from_o_ring(.443,70)
    force_injector_o_rings = sum_forces({force_outer_o_ring,force_film_o_ring,force_PT_o_ring}, 1)
    print(f"force_injector_o_rings: {force_injector_o_rings:.2f}")
    
    net_force_injector = sum_forces({force_manifold_pressure,force_injector_o_rings}, safety_factor)
    print(f"net_force_injector: {net_force_injector:.2f}")


    force_upper_o_ring = tensile_from_o_ring(5.44302,70)
    force_lower_o_ring = tensile_from_o_ring(5.44302,70)
    force_chamber_o_rings = sum_forces({force_upper_o_ring,force_lower_o_ring}, 1)
    print(f"force_chamber_o_rings: {force_chamber_o_rings:.2f}")

    net_force_chamber = sum_forces({force_manifold_pressure,force_chamber_o_rings}, safety_factor)
    print(f"net_force_chamber: {net_force_chamber:.2f}")


    force_pintle_o_ring = tensile_from_o_ring(1.505,70)
    print(f"force_pintle_o_ring: {force_pintle_o_ring:.2f}")

    net_force_pintle = sum_forces({force_pintle_pressure,force_pintle_o_ring}, safety_factor)
    print(f"net_force_pintle: {net_force_pintle:.2f}")

    flange_diameter = diameter(ED_ratio, bolt_major_diameter_main, nut_diameter_main, chamber_diameter, chamber_wall_thickness)
    print(f"flange_diameter: {flange_diameter:.2f}")


    lower_bound_preload_main = lower_preload(bolt_minor_diameter_main,proof_stress,preload_percent,torque_variation)
    upper_bound_preload_main = upper_preload(bolt_minor_diameter_main,proof_stress,preload_percent)
    print(f"lower_bound_preload_main: {lower_bound_preload_main:.2f}")
    print(f"upper_bound_preload_main: {upper_bound_preload_main:.2f}")

    lower_bound_preload_pintle = lower_preload(bolt_minor_diameter_pintle,proof_stress,preload_percent,torque_variation)
    upper_bound_preload_pintle = upper_preload(bolt_minor_diameter_pintle,proof_stress,preload_percent)
    print(f"lower_bound_preload_pintle: {lower_bound_preload_pintle:.2f}")
    print(f"upper_bound_preload_pintle: {upper_bound_preload_pintle:.2f}")


    calculated_number_of_bolts_injector = bolts(net_force_injector,lower_bound_preload_main)
    print(f"calculated_number_of_bolts_injector: {calculated_number_of_bolts_injector:.2f}")

    actual_number_bolts_injector = 6
    MOS_injector = MOS(calculated_number_of_bolts_injector, actual_number_bolts_injector)
    print(f"actual_number_bolts_injector: {actual_number_bolts_injector:.2f}")
    print(f"MOS_injector: {MOS_injector:.2f}")


    calculated_number_of_bolts_chamber= bolts(net_force_chamber,lower_bound_preload_main)
    print(f"calculated_number_of_bolts_chamber: {calculated_number_of_bolts_chamber:.2f}")

    actual_number_bolts_chamber = 6
    MOS_chamber = MOS(calculated_number_of_bolts_chamber, actual_number_bolts_chamber)
    print(f"actual_number_bolts_chamber: {actual_number_bolts_chamber:.2f}")
    print(f"MOS_chamber: {MOS_chamber:.2f}")


    calculated_number_of_bolts_pintle = bolts(net_force_pintle,lower_bound_preload_pintle)
    print(f"calculated_number_of_bolts_pintle: {calculated_number_of_bolts_pintle:.2f}")

    actual_number_bolts_pintle = 6
    MOS_pintle = MOS(calculated_number_of_bolts_pintle, actual_number_bolts_pintle)
    print(f"actual_number_bolts_pintle: {actual_number_bolts_pintle:.2f}")
    print(f"MOS_pintle: {MOS_pintle:.2f}")

