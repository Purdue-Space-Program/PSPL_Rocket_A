from math import exp, erfc

def heat_transfer_coefficient():
    Dt = 1 #diameter of engine throat
    mu = 1 #dynamic viscosity of the combustion gas
    Cp = 1 #specific heat of the combustion gas at constant pressure
    P0 = 1 #stagnation/chamber pressure (I think)
    g = 9.81 #gravitational constant
    Rt = 1 #radius of the throat curve
    Area_ratio = 1 #ratio of throat area to the local area at the point of interest
    gamma = 1 #specific heat ratio
    c_star = 1 #characteristic exhaust velocity
    Pr = 1 #Prandtl number of the combustion gas
    Twg = 1 #wall temperature (K)
    T0 = 1 #stagnation temperature (K)
    M = 1 #Mach number

    #The sigma term of the Bartz equation split into different terms
    sigma_parentheses1 = ((0.5 * (Twg / T0) * (1 + ((gamma - 1)/2) * (M**2))) + 0.5) ** 0.68
    sigma_parentheses2 = (1 + (((gamma - 1)/2) * (M**2))) ** 0.12
    sigma = (sigma_parentheses1 * sigma_parentheses2) ** (-1)

    #Bartz equations split into different terms
    heat_transfer_term1 = 0.026 / (Dt ** 0.2)
    heat_transfer_term2 = ((mu ** 0.2) * Cp) / (Pr ** 0.6)
    heat_transfer_term3 = ((P0 * g) / (c_star)) ** 0.8
    heat_transfer_term4 = (Dt / Rt) ** 0.1

    bartz_equation = heat_transfer_term1 * heat_transfer_term2 * heat_transfer_term3 * heat_transfer_term4 * (Area_ratio ** 0.9) * sigma

    return bartz_equation

def temperature_surface_calculation(heat_transfer_coefficient):

    alpha = 1 #thermal diffusivity of the chamber wall material
    T_infinity = 1 #combustion temperature (Kelvin)
    T_initial = 293 #room temperature in Kelvin
    t = 10 #burn time in seconds
    k = 1 #thermal conductivity of the chamber wall material

    #Surface temperature equation split into different terms
    Surface_temp_term1 = (-1) * exp(((heat_transfer_coefficient ** 2) * alpha * t) / (k**2))
    Surface_temp_term2 = erfc((heat_transfer_coefficient * ((alpha * t) ** 0.5)) / k)
    Surface_temp_term3 = T_infinity - T_initial 

    surface_temp = (Surface_temp_term1 * Surface_temp_term2 * Surface_temp_term3) + T_initial

    return surface_temp


def main():
    heat_transfer_coefficient = heat_transfer_coefficient()
    print("Heat flux onto the chamber wall:", heat_transfer_coefficient)

    temperature_surface = temperature_surface_calculation(heat_transfer_coefficient)
    print("Surface temperature of the chamber wall:", temperature_surface)









    #now that we have the heat transfer, we can calculate the heat flux 
    #I think this is done by dividing bartz equation by the fluid's specific heat multiplied by the temperature difference between the wall and the gas
'''

    gamma = 1 #specific heat of the fluid
    Tw = 1 #temperature of the wall
    Tg = 1 #temperature of the gas

    heat_flux =  heat_transfer_coefficient / (gamma * (Tw - Tg))
    print("Heat flux onto the chamber wall", heat_flux)

    surface_area = 1 #surface area of the chamber wall
    final_heatTransfer = heat_flux * surface_area
    print("final heat transfer onto the chamber wall", final_heatTransfer)

    #now to find the maximum wall temperature
    burn_time = 10 #burn time (seconds)
    wall_density = 1 #kg/(m^3)
    wall_specific_heat = 1 #specific heat of the wall material (J/(kg*K))
    wall_thickness = 1 #thickness of the wall (m)

    wall_mass = surface_area * wall_thickness * wall_density
    delta_T = final_heatTransfer * burn_time / (wall_mass * wall_specific_heat)

    '''




