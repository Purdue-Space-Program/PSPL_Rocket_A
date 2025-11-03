from math import exp, erfc, erf, pi
from scipy.optimize import fsolve

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
chamber_contour_csv_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "ChamberContour", "chamber_contour_meters"))

import constants as c

os.environ["CEA_USE_LEGACY"] = "1" # https://github.com/civilwargeeky/CEA_Wrap/issues/8
import CEA_Wrap as CEA

def main():
    M2IN = 39.3701
    IN2M = 1 / 39.3701


    #cylinder part of the chamber geometry parameters
    chamber_length = 11.167 * IN2M #chamber length (m)
    D_star = 2.3094013 * IN2M #throat diameter (m) # UPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATE
    A_star = pi * (D_star / 2)**2 #throat area (m^2)
    chamber_diameter = 6 * IN2M #chamber diameter (m)

    chamber_contour = np.loadtxt(chamber_contour_csv_path, delimiter=',')
    station_depths = chamber_contour[:, 0]
    station_inner_radii = chamber_contour[:, 1]

    station_areas = np.pi * (station_inner_radii**2)
    station_area_ratios = station_areas / A_star

    #initializing arrays to store Mach number, heat transfer coefficient, and surface temperature values for full chamber + nozzle
    Mach_total = np.zeros_like(station_area_ratios) #mach number at each axial position
    h_total = np.zeros_like(station_area_ratios) #heat transfer coefficient at each axial position
    Temp_surface_total = np.zeros_like(station_area_ratios) #surface temperature at each axial position

    initial_guess = 0.2

    #now calculating Mach number, heat transfer coefficient, and surface temperature at each position along the chamber length
    for station_index, A_ratio in enumerate(station_area_ratios):
        
        cea_results = RunCEA(150, "ethanol", "liquid oxygen", 1.0, sub = A_ratio[station_index])

        if station_index > 0:
            initial_guess = Mach_total[station_index - 1]
        else:
            initial_guess = 0.5
        
        M_local = calculating_MachNumber(gamma = cea_results["gamma"], area_ratio_value = A_ratio, initial_guess = initial_guess)
        Mach_total[station_index] = M_local

        #updating initial guess for next iteration
        initial_guess = M_local

        h_local = heat_transfer_coefficient(
            Dt = 2 * station_inner_radii[station_index],  # local diameter
            Rt = ((1.5 * 1.15 * IN2M) + (0.382 * 1.15 * IN2M)) / 2,     #radius of throat curve (m)
            Pr = cea_results["pran"], #Prandtl number of the combustion gas (n/a)
            gamma = cea_results["gamma"], #specific heat ratio of the combustion gas (n/a)
            c_star = cea_results["c_star"], #characteristic exhaust velocity (m/s)
            T0 = cea_results["t"], #stagnation temperature of the combustion gas ((K))
            Twg = recovery_temperature( #recovery temperature at the wall
                T_c = cea_results["t"],
                gamma = cea_results["gamma"],
                M = M_local,
                Pr = cea_results["pran"]
            ),
            Cp = cea_results["cp"] * 1000, #specific heat at constant pressure of the combustion gas (J/(kg*K))
            P0 = cea_results["p"] * 1e5, #chamber pressure (Pascals)
            mu = cea_results["visc"], #dynamic viscosity of the combustion gas (Pascal - seconds)
            M = Mach_total[station_index], #Mach number at the local axial point (no units)
            local_Area_ratio = A_ratio #area ratio at the local axial point (no units)

            )
        h_total[station_index] = h_local

        '''
        Temp_surface_total[station_index] = temperature_surface_calculation(
            heat_transfer_coefficient_value = h_total[station_index],
            T_infinity = cea_results["c_t"], #chamber temperature (K)
            k = cea_results["c_cond"] #conductivity of the combustion gas in the chamber (W/(m*K))
        )
        '''
        Temp_surface_total[station_index] = temperature_surface_calculation(
            heat_transfer_coefficient_value = h_total[station_index],
            axial_position = station_depths[station_index],
            T_infinity = cea_results["t"], #chamber temperature (K)
            k = 50 #thermal conductivity of the chamber wall material (W/(m*K))
        )
            
    #printing results

    #printing axial positions vs surface temp plot
    plt.figure()
    plt.plot(station_depths * M2IN, Temp_surface_total)
    plt.xlabel("Axial Position Relative to Throat (in) ")
    plt.ylabel("Surface Temperature (K) ")
    plt.title("Surface temperature vs Axial Position")
    plt.grid(True)
    plt.show()

    #printing axial position vs heat transfer coefficient
    plt.figure()
    plt.plot(station_depths * M2IN, h_total)
    plt.xlabel("Axial Position Relative to Throat (in) ")
    plt.ylabel("Heat Transfer Coefficient (W/m^2 K) ")
    plt.title("Heat Transfer Coefficient vs Axial Position")
    plt.grid(True)
    plt.show()

    print("Maximum Surface Temperature (K): ", max(Temp_surface_total))
    print("Maximum Heat Transfer Coefficient (W/m^2 K): ", max(h_total))
    

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

    exit_pressure = 15 # [psi]
    pressure_ratio = chamber_pressure / exit_pressure # assume exit pressure is a constantly at the pressure of air a bit above sea level

    rocket = CEA.RocketProblem(
        pressure =       chamber_pressure,
        pip =            pressure_ratio, # pip is "Pressure ratio of chamber pressure to exit pressure." github.com/civilwargeeky/CEA_Wrap/blob/main/README.md#rocket-problem-constructor-additional-parameters
        materials =      [CEA_fuel_name, CEA_oxidizer_name],
        o_f =            OF_Ratio,
        pressure_units = "psi",
        analysis_type= "frozen"
    )

    cea_results = rocket.run()
    return{
        "p": cea_results.p, #chamber pressure (Bar)
        "c_star": cea_results.cstar, #characteristic exhaust velocity in (m/s)
        "pran": cea_results.pran, #Prandtl number of combustion gas at chamber (no units)
        "gamma": cea_results.gamma, #specific heat ratio of combustion gas at chamber (no units)
        "t": cea_results.t, #stagnation temperature at chamber (K)
        "cp": cea_results.cp, #specific heat at constant pressure of combustion gas (kJ/kg*K)
        "visc": cea_results.visc, #dynamic viscosity of combustion gas in the chamber (Pascal - seconds)
        "cond": cea_results.cond, #conductivity of combustion gas in the chamber (W/m*K)
    }
    
def recovery_temperature(T_c, gamma, M, Pr):

    T_r = T_c * (1 + (((Pr*(gamma - 1)) / 2) * M**2))

    return T_r


def calculating_MachNumber(gamma, area_ratio_value, initial_guess):

    def f(M):
        Mach_function_part1 = ((gamma + 1)/2)**(-(gamma + 1)/(2*(gamma-1)))
        Mach_function_part2 = (1/M) * ((1 + (((gamma - 1)/2) * M**2)) ** ((gamma + 1)/(2*(gamma-1))))
        Mach_function = (Mach_function_part1 * Mach_function_part2) - area_ratio_value
        return abs(Mach_function)
    
    
    M_solution = fsolve(f, initial_guess)
    return float(M_solution)

'''
def heat_transfer_coefficient(Dt, Rt, Pr, gamma, c_star, T0, Twg, Cp, P0, mu, M, local_Area_ratio):

    #The sigma term of the Bartz equation split into different terms
    sigma_parentheses1 = ((0.5 * (Twg / T0) * (1 + (((gamma - 1) * M**2)/2))) + 0.5) ** 0.68
    sigma_parentheses2 = (1 + (((gamma - 1)/2) * (M**2))) ** 0.12
    sigma = 1 / (sigma_parentheses1 * sigma_parentheses2)

    #Bartz equations split into different terms
    heat_transfer_term1 = 0.026 / (Dt ** 0.2)
    heat_transfer_term2 = ((mu ** 0.2) * Cp) / (Pr ** 0.6)
    heat_transfer_term3 = (P0 / (c_star)) ** 0.8
    heat_transfer_term4 = (Dt / Rt) ** 0.1

    bartz_equation = heat_transfer_term1 * heat_transfer_term2 * heat_transfer_term3 * heat_transfer_term4 * (local_Area_ratio ** 0.9) * sigma

    return bartz_equation


def temperature_surface_calculation(heat_transfer_coefficient_value, T_infinity, k):
    
    alpha = 3.6e-4 #thermal diffusivity of the chamber wall material ((m^2)/s)
    T_initial = 293 #room temperature (K)
    t = 20 #burn time (sec)

    #Surface temperature equation split into different terms
    Surface_temp_term0 = (heat_transfer_coefficient_value ** 2) * alpha * t
    Surface_temp_term1 = exp(-(Surface_temp_term0) / (k**2))
    Surface_temp_term2 = erfc((heat_transfer_coefficient_value * ((alpha * t) ** 0.5)) / k)
    Surface_temp_term3 = T_infinity - T_initial 
    surface_temp = (Surface_temp_term1 * Surface_temp_term2 * Surface_temp_term3) + T_initial 

    return surface_temp        

'''

def heat_transfer_coefficient(Dt, Rt, Pr, gamma, c_star, T0, Twg, Cp, P0, mu, M, local_Area_ratio):
    # sigma (Bartz correction)
    sigma_parentheses1 = ((0.5 * (Twg / T0) * (1 + (((gamma - 1) * M**2)/2))) + 0.5) ** 0.68
    sigma_parentheses2 = (1 + (((gamma - 1)/2) * (M**2))) ** 0.12
    sigma = 1.0 / (sigma_parentheses1 * sigma_parentheses2)

    # Bartz core terms 
    heat_transfer_term1 = 0.026 / (Dt ** 0.2)
    heat_transfer_term2 = ((mu ** 0.2) * Cp) / (Pr ** 0.6)
    heat_transfer_term3 = (P0 / (c_star)) ** 0.8
    heat_transfer_term4 = (Dt / Rt) ** 0.1

    # using inverse of local_Area_ratio so heat goes up as A goes down 
    area_factor = (1.0 / local_Area_ratio) ** 0.9

    bartz_equation = heat_transfer_term1 * heat_transfer_term2 * heat_transfer_term3 * heat_transfer_term4 * area_factor * sigma

    return bartz_equation

def temperature_surface_calculation(heat_transfer_coefficient_value, axial_position, T_infinity, k = 167):
    Ti = 294 #K, initial temperature of the chamber wall
    alpha = 1.5e-5 #thermal diffusivity of the chamber wall material 
    t = 2.24 #sec, burn time

    term_conduction = k / ((pi * alpha * t) ** 0.5)     # conduction resistance term
    Ts = (heat_transfer_coefficient_value * T_infinity + term_conduction * Ti) / (heat_transfer_coefficient_value + term_conduction)

    return Ts


if __name__ == "__main__":
    main()
