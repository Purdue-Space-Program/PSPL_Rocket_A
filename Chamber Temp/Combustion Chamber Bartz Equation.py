from math import exp, erfc
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import coding_utils.constants as c
import vehicle_scripts.numpy_ndarray_handler as NNH

os.environ["CEA_USE_LEGACY"] = "1" # https://github.com/civilwargeeky/CEA_Wrap/issues/8
import CEA_Wrap as CEA

def heat_transfer_coefficient():
    Dt = 0.0498 #diameter of engine throat
    mu = 4e-5 #dynamic viscosity of the combustion gas (Pa*s)
    Cp = 2500 #specific heat of the combustion gas at constant pressure(J/kg*K)
    P0 = 3e6 #stagnation/chamber pressure (I think) (Pa)
    g = 9.81 #gravitational constant
    Rt = 0.05 #radius of the throat curve
    Area_ratio = 1 #ratio of throat area to the local area at the point of interest
    gamma = 1.2 #specific heat ratio, typical value for hot combustion gases
    c_star = 1248 #characteristic exhaust velocity
    Pr = 0.7 #Prandtl number of the combustion gas
    Twg = 800 #wall temperature (K) because steel can withstand up to 1100 K but safety margin
    T0 = 3500 #stagnation temperature (K)
    M = 1 #Mach number at the throat

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

def temperature_surface_calculation(heat_transfer_coefficient_value):

    alpha = 4e-6 #thermal diffusivity of the chamber wall material
    T_infinity = 3500 #combustion temperature (Kelvin)
    T_initial = 293 #room temperature in Kelvin
    t = 10 #burn time in seconds
    k = 16 #thermal conductivity of the chamber wall material (W/mK)

    #Surface temperature equation split into different terms
    Surface_temp_term0 = (heat_transfer_coefficient_value ** 2) * alpha * t
    Surface_temp_term1 = (-1) * exp((Surface_temp_term0) / (k**2))
    Surface_temp_term2 = erfc((heat_transfer_coefficient_value * ((alpha * t) ** 0.5)) / k)
    Surface_temp_term3 = T_infinity - T_initial 

    surface_temp = (Surface_temp_term1 * Surface_temp_term2 * Surface_temp_term3) + T_initial

    return surface_temp

'''
def minimum_wall_thickness():
    #using the formula (PD)/(2(SE+PY))
    P = 0.000435113 #chamber pressure (psi)
    D = 0.0498 #outer diameter of engine throat
    S = #allowed stress value which depends on the specific grade of stainless steel and design temp
    E = 1 #longitudinal weld join quality factor, depdends if material has been radiographed
    Y = 0.7 #wall thickness coefficient, stainless steel for above 900F

    thickness = (P * D) / (2 * ((S * E) + (P * Y)))
    return thickness
'''

def main():
    heat_transfer_coefficient_value = heat_transfer_coefficient()
    print("Heat flux (h) onto the chamber wall:", heat_transfer_coefficient_value)

    temperature_surface = temperature_surface_calculation(heat_transfer_coefficient_value)
    print("Surface temperature of the chamber wall (in Kelvin):", temperature_surface)

    '''
    thickness = minimum_wall_thickness()
    print("Minimum wall thickness of the chamber (in units idk yet):", thickness)
    '''

if __name__ == "__main__":
    main()




