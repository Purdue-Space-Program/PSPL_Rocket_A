from math import exp, erfc
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

os.environ["CEA_USE_LEGACY"] = "1" # https://github.com/civilwargeeky/CEA_Wrap/issues/8
import CEA_Wrap as CEA

def RunCEA(
    chamber_pressure,
    fuel_name,
    oxidizer_name,
    OF_Ratio,
):
    # convert regular string for propellants to what CEA_wrap uses
    if fuel_name == "ethanol":
        CEA_fuel_name = CEA.Fuel("C2H5OH(L)", temp=290)
    elif fuel_name == "kerosene":
        CEA_fuel_name = CEA.Fuel("Jet-A(L)", temp=290)
    elif fuel_name == "ipa":
        CEA_fuel_name = CEA.Fuel("C3H8O,2propanol", temp=290)
    else:
        raise ValueError(f"{fuel_name} not supported")

    if oxidizer_name == "liquid oxygen":
        CEA_oxidizer_name = CEA.Oxidizer("O2(L)", temp=90) # 90 K is temperature of oxidizer upon injection into combustion (same as copperhead's sizing)
    else:
        raise ValueError(f"{oxidizer_name} not supported")


    pressure_ratio = chamber_pressure / (15 * 6894.76) # assume exit pressure is a constantly at the pressure of air a bit above sea level

    rocket = CEA.RocketProblem(
        pressure =       chamber_pressure * (1 / 6894.76),
        pip =            pressure_ratio, # pip is "Pressure ratio of chamber pressure to exit pressure." github.com/civilwargeeky/CEA_Wrap/blob/main/README.md#rocket-problem-constructor-additional-parameters
        materials =      [CEA_fuel_name, CEA_oxidizer_name],
        o_f =            OF_Ratio,
        pressure_units = "psi",
    )

    cea_results = rocket.run()
    return{
        "c_p": cea_results.c_p, #chamber pressure
        "c_star": cea_results.cstar, #characteristic exhaust velocity
        "c_pran": cea_results.c_pran, #Prandtl number of combustion gas
        "gamma": cea_results.gamma, #specific heat ratio of combustion gas
        "c_t": cea_results.c_t #stagnation temperature
    }

def heat_transfer_coefficient(Pr, gamma, c_star, T0):
    Dt = 0.0498 #diameter of engine throat
    mu = 4e-5 #dynamic viscosity of the combustion gas (Pa*s)
    Cp = 2.2996 #specific heat of the combustion gas at constant pressure(J/kg*K)
    P0 = 3e6 #stagnation/chamber pressure (I think) (Pa)
    g = 9.81 #gravitational constant
    Rt = 0.05 #radius of the throat curve
    Area_ratio = 3 #ratio of throat area to the local area at the point of interest (contraction ratio)
    Twg = 800 #wall temperature (K) because steel can withstand up to 1100 K but safety margin
    M = 2.16 #Mach number at the throat

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

def temperature_surface_calculation(heat_transfer_coefficient_value, T_infinity):

    alpha = 3.6e-6 #thermal diffusivity of the chamber wall material
    T_initial = 293 #room temperature in Kelvin
    t = 2.24 #burn time in seconds
    k = 16.3 #thermal conductivity of the chamber wall material (W/mK)

    #Surface temperature equation split into different terms
    Surface_temp_term0 = (heat_transfer_coefficient_value ** 2) * alpha * t
    Surface_temp_term1 = exp(-(Surface_temp_term0) / (k**2))
    Surface_temp_term2 = erfc((heat_transfer_coefficient_value * ((alpha * t) ** 0.5)) / k)
    Surface_temp_term3 = T_infinity - T_initial 

    surface_temp = (Surface_temp_term1 * Surface_temp_term2 * Surface_temp_term3) + T_initial 

    return surface_temp

'''
def minimum_wall_thickness():
    #using the formula (PD)/(2(SE+PY))
    P = 150 #chamber pressure (psi)
    D = 0.0498 #outer diameter of engine throat
    S = #allowed stress value which depends on the specific grade of stainless steel and design temp
    E = 1 #longitudinal weld join quality factor, depdends if material has been radiographed
    Y = 0.7 #wall thickness coefficient, stainless steel for above 900F

    thickness = (P * D) / (2 * ((S * E) + (P * Y)))
    return thickness
'''

def main():
    PSI2PA = 6894.76 # [Pa/psi] Conversion factor from psi to Pa

    cea_results = RunCEA(150 * PSI2PA, "ethanol", "liquid oxygen", 1.0)

    heat_transfer_coefficient_value = heat_transfer_coefficient(
        Pr = cea_results["c_pran"],
        gamma = cea_results["gamma"],
        c_star = cea_results["c_star"],
        T0 = cea_results["c_t"],
    )
    
    print("Heat flux (h) onto the chamber wall:", heat_transfer_coefficient_value)

    temperature_surface = temperature_surface_calculation(
        heat_transfer_coefficient_value,
        T_infinity = cea_results["c_t"])
    print("Surface temperature of the chamber wall (in Kelvin):", temperature_surface)

    '''
    thickness = minimum_wall_thickness()
    print("Minimum wall thickness of the chamber (in units idk yet):", thickness)
    '''



if __name__ == "__main__":
    main()





