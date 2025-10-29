from math import exp, erfc, pi
from scipy.optimize import fsolve
import sys
import os
import numpy as np
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
        "c_p": cea_results.c_p, #chamber pressure (Bar)
        "c_star": cea_results.cstar, #characteristic exhaust velocity (m/s)
        "c_pran": cea_results.c_pran, #Prandtl number of combustion gas (no units)
        "gamma": cea_results.gamma, #specific heat ratio of combustion gas (no units)
        "c_t": cea_results.c_t, #stagnation temperature (K)
        "c_cp": cea_results.c_cp, #specific heat at constant pressure of combustion gas (kJ/kg*K)
        "c_visc": cea_results.c_visc, #dynamic viscosity of combustion gas (Pa*s)
        "c_cond": cea_results.c_cond, #conductivity of combustion gas in the chamber (W/m*K)
        "mach": cea_results.mach, #Mach number at the nozzle exit (no units)
    }

def calculating_MachNumber(gamma, area_ratio_value, initial_guess = 0.2):

    def f(M):
        Mach_function_part1 = ((gamma + 1)/2)**(-(gamma + 1)/(2*(gamma-1)))
        Mach_function_part2 = (1/M) * ((1 + (((gamma - 1)/2) * M**2)) ** ((gamma + 1)/(2*(gamma-1))))
        Mach_function = (Mach_function_part1 * Mach_function_part2) - area_ratio_value
        return Mach_function
    
    
    M_solution = fsolve(f, initial_guess)
    return M_solution


def heat_transfer_coefficient(Pr, gamma, c_star, T0, Cp, P0, mu, M, local_Area_ratio):
    Dt = 0.1524 #diameter of chamber (m)
    Rt = 0.05 #radius of the throat curve 
    Twg = 800 #wall temperature (K) because steel can withstand up to 1100 K but safety margin

    #The sigma term of the Bartz equation split into different terms
    sigma_parentheses1 = ((0.5 * (Twg / T0) * (1 + ((gamma - 1)/2) * (M**2))) + 0.5) ** 0.68
    sigma_parentheses2 = (1 + (((gamma - 1)/2) * (M**2))) ** 0.12
    sigma = (sigma_parentheses1 * sigma_parentheses2) ** (-1)

    #Bartz equations split into different terms
    heat_transfer_term1 = 0.026 / (Dt ** 0.2)
    heat_transfer_term2 = ((mu ** 0.2) * Cp) / (Pr ** 0.6)
    heat_transfer_term3 = (P0 / (c_star)) ** 0.8
    heat_transfer_term4 = (Dt / Rt) ** 0.1

    bartz_equation = heat_transfer_term1 * heat_transfer_term2 * heat_transfer_term3 * heat_transfer_term4 * (local_Area_ratio ** 0.9) * sigma

    return bartz_equation

def temperature_surface_calculation(heat_transfer_coefficient_value, T_infinity, k):
    
    alpha = 3.6e-6 #thermal diffusivity of the chamber wall material ((m^2)/s)
    T_initial = 293 #room temperature (K)
    t = 2.24 #burn time (sec)

    #Surface temperature equation split into different terms
    Surface_temp_term0 = (heat_transfer_coefficient_value ** 2) * alpha * t
    Surface_temp_term1 = exp(-(Surface_temp_term0) / (k**2))
    Surface_temp_term2 = erfc((heat_transfer_coefficient_value * ((alpha * t) ** 0.5)) / k)
    Surface_temp_term3 = T_infinity - T_initial 
    surface_temp = (Surface_temp_term1 * Surface_temp_term2 * Surface_temp_term3) + T_initial 

    return surface_temp


def main():
    PSI2PA = 6894.76 # [Pa/psi] Conversion factor from psi to Pa

    cea_results = RunCEA(150 * PSI2PA, "ethanol", "liquid oxygen", 1.0)

    #Chamber geometery parameters
    chamber_length = 20 #chamber length (m)
    dx = 0.001 #increments of 1mm
    D_star = 1 #throat diamater (m)
    Astar = pi * (D_star / 2)**2 #throat area (m^2)
    chamber_diameter = 3.5 #chamber diameter (m)
    chamber_area = pi * (chamber_diameter/2)**2 #chamber area (m^2)

    #linearly interpolating area along chamber length (hopefully it's a straight line)
    x_positions = np.arange(0,chamber_length + dx, dx) #position along the chamber length (m)
    area_values = np.linspace(Astar, chamber_area, len(x_positions)) #local areas (m^2)
    area_ratios = area_values / Astar #area ratio A/A* (no units)

    #initializing arrays to store Mach number, heat transfer coefficient, and surface temperature values
    Mach_array = np.zeros_like(area_ratios) #mach number at each axial position
    h_array = np.zeros_like(area_ratios) #heat transfer coefficient at each axial position
    Temp_surface_array = np.zeros_like(area_ratios) #surface temperature at each axial position

    initial_guess = 2

    #now calculating Mach number, heat transfer coefficient, and surface temperature at each position along the chamber length
    for i, A_ratio in enumerate(area_ratios):
        
        M_local = calculating_MachNumber(gamma = cea_results["gamma"], area_ratio_value = A_ratio, initial_guess = initial_guess)
        Mach_array[i] = M_local

        #updating initial guess for next iteration
        initial_guess = M_local

        h_local = heat_transfer_coefficient(
            Pr = cea_results["c_pran"], #Prandlt number of the combustion gas       
            gamma = cea_results["gamma"], #specific heat ratio of the combustion gas
            c_star = cea_results["c_star"], #characteristic exhaust velocity (m/s)
            T0 = cea_results["c_t"], #stagnation temperature of the combustion gas ((K))
            Cp = cea_results["c_cp"], #specific heat at constant pressure of the combustion gas
            P0 = cea_results["c_p"] * 100000, #chamber pressure (Bar)
            mu = cea_results["c_visc"], #dynamic viscosity of the combustion gas (Pa
            M = Mach_array[i], #Mach number at the local axial point (no units)
            local_Area_ratio = area_ratios[i] #area ratio at the local axial point (no units)

        )
        h_array[i] = h_local

        Temp_surface_array[i] = temperature_surface_calculation(
            heat_transfer_coefficient_value = h_array[i],
            T_infinity = cea_results["c_t"], #chamber temperature (K)
            k = cea_results["c_cond"] #conductivity of the combustion gas in the chamber (W/(m*K))
        )


if __name__ == "__main__":
    main()





