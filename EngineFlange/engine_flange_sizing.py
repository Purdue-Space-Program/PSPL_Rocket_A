pi = 3.14159

def tensile_from_chamber(diameter, pressure):
    return (diameter*diameter/4.0*pi) * pressure

def tensile_from_o_ring(diameter, compression):
    return (diameter*pi) * compression

def proof(safety, strength):
    return safety*strength

def lower_preload(diameter, proof):
    area = pi * (diameter*diameter/4)
    return (((area*proof)*0.75)/1.25)*0.75

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

def main ():
    proof_stress_safety_factor = 0.8
    steel_tensile_strength = 42100
    proof_stress = proof(proof_stress_safety_factor,steel_tensile_strength)
    safety_factor = 1.4
    #Pathfinder
    bolt_diameter_plates = 0.1875
    ED_ratio = 1.5
    nut_diameter_plates = 0.375
    bolt_diameter_pintle = 0.125
    chamber_diameter = 4.9
    chamber_wall_thickness = 0.25
    chamber_pressure = 250
    #CMS
    #bolt_diameter_plates = 0.25
    #chamber_diameter = 5.05
    #chamber_pressure = 200

    #Pathfinder
    tensile_force_from_chamber_plates = tensile_from_chamber(chamber_diameter,chamber_pressure)
    tensile_force_from_outer_o_ring = tensile_from_o_ring(5.19301,70)
    tensile_force_from_chamber_o_ring = tensile_from_o_ring(5.19302,70)
    tensile_force_from_film_o_ring = 0 #tensile_from_o_ring(3.984,70)
    tensile_force_from_manifold_o_ring = 0 #tensile_from_o_ring(3.512,70)

    #CMS
    #tensile_force_from_chamber_plates = tensile_from_chamber(chamber_diameter,chamber_pressure)
    #tensile_force_from_outer_o_ring = 0
    #tensile_force_from_chamber_o_ring = tensile_from_o_ring(5.19301,70)
    #tensile_force_from_film_o_ring = 0 
    #tensile_force_from_manifold_o_ring = 0 

    forces_plates = {tensile_force_from_chamber_plates,tensile_force_from_outer_o_ring,tensile_force_from_chamber_o_ring,tensile_force_from_film_o_ring,tensile_force_from_manifold_o_ring}
    net_force_plates = sum_forces(forces_plates, safety_factor)
    lower_bound_preload_plates = lower_preload(bolt_diameter_plates,proof_stress)
    number_of_bolts_plates = bolts(net_force_plates,lower_bound_preload_plates)
    print(number_of_bolts_plates)

    tensile_force_from_chamber_pintle = tensile_from_chamber(.98,500)
    tensile_force_from_pintle_o_ring = tensile_from_o_ring(1.25,70)

    flange_diameter = diameter(ED_ratio, bolt_diameter_plates, nut_diameter_plates, chamber_diameter, chamber_wall_thickness)
    print(flange_diameter)

    forces_pintle = {tensile_force_from_chamber_pintle,tensile_force_from_pintle_o_ring}
    net_force_pintle = sum_forces(forces_pintle, safety_factor)
    lower_bound_preload_pintle = lower_preload(bolt_diameter_pintle,proof_stress)
    number_of_bolts_pintle = bolts(net_force_pintle,lower_bound_preload_pintle)
    print(number_of_bolts_pintle)

main()
