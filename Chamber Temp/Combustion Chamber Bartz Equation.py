from math import exp, erfc, erf, pi
from scipy.optimize import fsolve
from dataclasses import dataclass

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
chamber_contour_csv_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "ChamberContour", "chamber_contour_meters.csv"))


import vehicle_parameters as vp
import constants as c

os.environ["CEA_USE_LEGACY"] = "1" # https://github.com/civilwargeeky/CEA_Wrap/issues/8
import CEA_Wrap as CEA


def main():

    #cylinder part of the chamber geometry parameters
    D_star = vp.parameters.chamber_throat_diameter #throat diameter (m) 
    A_star = pi * (D_star / 2)**2 #throat area (m^2)

    chamber_contour = np.loadtxt(chamber_contour_csv_path, delimiter=',')
    station_depths = chamber_contour[:, 0]
    station_inner_radii = chamber_contour[:, 2]

    #aligning axial positions to the throat
    throat_index = np.argmin(station_inner_radii)
    x_relative = station_depths - station_depths[throat_index]

    station_areas = np.pi * (station_inner_radii**2)
    station_area_ratios = station_areas / A_star

    #initializing arrays to store Mach number, heat transfer coefficient, and surface temperature values for full chamber + nozzle
    #Mach_total = np.zeros_like(station_area_ratios) #mach number at each axial position
    Mach_total = np.full_like(station_area_ratios, np.nan, dtype = float)
    #h_total = np.zeros_like(station_area_ratios) #heat transfer coefficient at each axial position
    h_total = np.full_like(station_area_ratios, np.nan, dtype = float)
    #Temp_surface_total = np.zeros_like(station_area_ratios) #surface temperature at each axial position
    Temp_surface_total = np.full_like(station_area_ratios, np.nan, dtype = float)


    cea_results = RunCEA(vp.parameters.chamber_pressure, "ethanol", "liquid oxygen", 1.0)
    gamma_loc = cea_results["gamma"]
    Pr_loc = cea_results["pran"]
    cp_loc = cea_results["cp"]
    #idk what units chamber pressure is in
    if cp_loc < 100:
        cp_loc *= 1000 #avoiding unrealistically low cp values from CEA
    visc_loc = cea_results["visc"]
    t_loc = cea_results["t"]
    P_loc = cea_results["p"] * 1e5
    if P_loc < 5000:
        p_loc *= 6894.76


    #now calculating Mach number, heat transfer coefficient, and surface temperature at each position along the chamber length
    for station_index, A_ratio in enumerate(station_area_ratios):

        #new stuff
        branch = "subsonic" if x_relative[station_index] <= 0 else "supersonic"
        initial_guess = Mach_total[station_index -1] if station_index > 0 else (0.5 if branch == "subsonic" else 2.0) 
        M_local = calculating_MachNumber(gamma = gamma_loc, area_ratio_value = A_ratio, initial_guess = initial_guess, branch = branch)
        Mach_total[station_index] = M_local 

        
        '''
        if station_index > 0:
            initial_guess = Mach_total[station_index - 1]
        else:
            initial_guess = 0.5
        M_local = calculating_MachNumber(gamma = gamma_loc, area_ratio_value = A_ratio, initial_guess = initial_guess)
        Mach_total[station_index] = M_local 
        '''
        Dt = 2 * station_inner_radii[station_index]  # local diameter

        if Dt == 0:
            Mach_total[station_index] = 0
            h_total[station_index] = 0
            Temp_surface_total[station_index] = 294 #K, initial temperature of the chamber wall
            continue
                
        #updating initial guess for next iteration
        initial_guess = M_local  

        h_local = heat_transfer_coefficient(
            Dt = Dt,  # local diameter
            Rt = ((1.5 * 1.15 * c.IN2M) + (0.382 * 1.15 * c.IN2M)) / 2,     #radius of throat curve (m)
            Pr = Pr_loc, #Prandtl number of the combustion gas (n/a)
            gamma = gamma_loc, #specific heat ratio of the combustion gas (n/a)
            c_star = cea_results["c_star"], #characteristic exhaust velocity (m/s)
            T0 = t_loc, #stagnation temperature of the combustion gas ((K))
            Twg = recovery_temperature( #recovery temperature at the wall
                T_c = t_loc,
                gamma = gamma_loc,
                M = M_local,
                Pr = Pr_loc
            ),
            Cp = cp_loc * 1000, #specific heat at constant pressure of the combustion
            P0 = P_loc, #chamber pressure (Pascals)
            mu = visc_loc, #dynamic viscosity of the combustion gas (Pascal - seconds)
            M = M_local, #Mach number at the local axial point (no units)
            local_Area_ratio = A_ratio #area ratio at the local axial point (no units)
            )
        h_total[station_index] = h_local

        Temp_surface_total[station_index] = temperature_surface_calculation(
            heat_transfer_coefficient_value = h_total[station_index],
            axial_position = station_depths[station_index],
            T_infinity = t_loc, #chamber temperature (K)
            k = 50, #thermal conductivity of the chamber wall material (W/(m*K))
            t = vp.parameters.burn_time #s, burn time
        )

        #plots
    plt.figure()
    plt.plot(station_depths * c.M2IN, Temp_surface_total)
    plt.xlabel("Axial Position Relative to Throat (in) ")
    plt.ylabel("Surface Temperature (K) ")
    plt.title("Surface temperature vs Axial Position")
    plt.grid(True)
    plt.show()

    plt.figure()
    plt.plot(station_depths * c.M2IN, h_total)
    plt.xlabel("Axial Position Relative to Throat (in) ")
    plt.ylabel("Heat Transfer Coefficient (W/m^2 K) ")
    plt.title("Heat Transfer Coefficient vs Axial Position")
    plt.grid(True)
    plt.show()
    

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
        analysis_type= "frozen",
    )

    cea_results = rocket.run()

    #print ("CEA result keys:", list(cea_results.keys()))

    return{
        "gamma": cea_results.gamma,
        "t":     cea_results.t,
        "visc":  cea_results.visc,
        "cp":    cea_results.cp,
        "pran":  cea_results.pran,
        "c_star": cea_results.cstar,
        "p":     cea_results.p,
    }
    
def recovery_temperature(T_c, gamma, M, Pr):

    T_r = T_c * (1 + (((Pr*(gamma - 1)) / 2) * M**2))

    return T_r


def calculating_MachNumber(gamma, area_ratio_value, initial_guess, branch):

    def f(M):
        Mach_function_part1 = ((gamma + 1)/2)**(-(gamma + 1)/(2*(gamma-1)))
        Mach_function_part2 = (1/M) * ((1 + (((gamma - 1)/2) * M**2)) ** ((gamma + 1)/(2*(gamma-1))))
        Mach_function = (Mach_function_part1 * Mach_function_part2) - area_ratio_value
        return Mach_function

    guess = initial_guess
    if branch == "subsonic:":
        guess = min(0.8, max(0.05, guess))

    else:
        guess = max(1.1, guess)

    M_solution = fsolve(f, initial_guess)
    return float(M_solution[0])


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
    area_factor = (1/local_Area_ratio) ** 0.9

    bartz_equation = heat_transfer_term1 * heat_transfer_term2 * heat_transfer_term3 * heat_transfer_term4 * area_factor * sigma

    return bartz_equation

def temperature_surface_calculation(heat_transfer_coefficient_value, axial_position, T_infinity, k = 167, t = vp.parameters.burn_time):
    Ti = 294 #K, initial temperature of the chamber wall
    alpha = 1.5e-5 #thermal diffusivity of the chamber wall material (m^2/s)

    term_conduction = k / ((pi * alpha * t) ** 0.5)     # conduction resistance term
    Ts = (heat_transfer_coefficient_value * T_infinity + term_conduction * Ti) / (heat_transfer_coefficient_value + term_conduction)

    return Ts


if __name__ == "__main__":
    main()
