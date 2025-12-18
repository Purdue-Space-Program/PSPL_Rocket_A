def tensile_from_chamber(diameter, pressure):
    return (diameter*diameter/4.0*3.14159) * pressure

def tensile_from_o_ring(diameter, compression):
    return (diameter*3.14159) * compression

def proof(safety, strength):
    return safety*strength

def lower_preload(diameter, proof):
    area = 3.14159 * (diameter*diameter/4)
    return (((area*proof)*0.75)/1.25)*0.75

def bolts(force, preload):
    return force/preload

def main ():
    proof_stress_safety_factor = 0.8
    steel_tensile_strength = 42100
    proof_stress = proof(proof_stress_safety_factor,steel_tensile_strength)

    #Pathfinder
    tensile_force_from_chamber_faceplate = tensile_from_chamber(4.9,500)
    tensile_force_from_outer_o_ring = tensile_from_o_ring(5.19301,70)
    tensile_force_from_chamber_o_ring = tensile_from_o_ring(5.19302,70)
    tensile_force_from_film_o_ring = 0 #tensile_from_o_ring(3.984,70)
    tensile_force_from_manifold_o_ring = 0 #tensile_from_o_ring(3.512,70)
    force_safety_factor_faceplate = 1.5

    #CMS
    #tensile_force_from_chamber_faceplate = tensile_from_chamber(5.05,200)
    #tensile_force_from_outer_o_ring = 0
    #tensile_force_from_chamber_o_ring = tensile_from_o_ring(5.19301,70)
    #tensile_force_from_film_o_ring = 0 
    #tensile_force_from_manifold_o_ring = 0 
    #force_safety_factor_faceplate = 1.5

    forces_faceplate = {tensile_force_from_chamber_faceplate,tensile_force_from_outer_o_ring,tensile_force_from_chamber_o_ring,tensile_force_from_film_o_ring,tensile_force_from_manifold_o_ring}
    net_force_faceplate = sum(forces_faceplate)*force_safety_factor_faceplate
    #print(net_force_faceplate)
    
    #Pathfinder
    bolt_diameter_faceplate = 0.25

    #CMS
    #bolt_diameter_faceplate = 0.25

    lower_bound_preload_faceplate=lower_preload(bolt_diameter_faceplate,proof_stress)
    #print(lower_bound_preload_faceplate)

    number_of_bolts_faceplate = bolts(net_force_faceplate,lower_bound_preload_faceplate)
    print(number_of_bolts_faceplate)

    tensile_force_from_chamber_pintle = tensile_from_chamber(.98,250)
    tensile_force_from_pintle_o_ring = tensile_from_o_ring(1.25,70)
    force_safety_factor_pintle = 1.5

    forces_pintle = {tensile_force_from_chamber_pintle,tensile_force_from_pintle_o_ring}
    net_force_pintle = sum(forces_pintle)*force_safety_factor_pintle
    #print(net_force_pintle)
    
    bolt_diameter_pintle = 0.125

    lower_bound_preload_pintle=lower_preload(bolt_diameter_pintle,proof_stress)
    #print(lower_bound_preload_pintle)

    number_of_bolts_pintle = bolts(net_force_pintle,lower_bound_preload_pintle)
    print(number_of_bolts_pintle)

main()
