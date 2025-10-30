from math import exp, erfc, pi
from scipy.optimize import fsolve
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c

chamber_contour_csv_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "ChamberContour", "chamber_contour_meters"))


os.environ["CEA_USE_LEGACY"] = "1" # https://github.com/civilwargeeky/CEA_Wrap/issues/8
import CEA_Wrap as CEA





def main():

    cea_results = RunCEA(150, "ethanol", "liquid oxygen", 1.0)

    #cylinder part of the chamber geometry parameters
    chamber_length = 11.167 * c.IN2M #chamber length (m)
    dx = 0.001 #increments of 1mm
    D_star = 2.3094013 * c.IN2M #throat diameter (m) # UPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATEUPDATE
    A_star = pi * (D_star / 2)**2 #throat area (m^2)
    chamber_diameter = 6 * c.IN2M #chamber diameter (m)
    chamber_area = pi * ((chamber_diameter/2)**2) #chamber area (m^2)


    chamber_contour = np.loadtxt(chamber_contour_csv_path, delimiter=',')
    station_depths = chamber_contour[:, 0]
    station_inner_radii = chamber_contour[:, 0]

    station_areas = np.pi * (station_inner_radii**2)
    station_area_ratios = station_areas / A_star

    





    # #linearly interpolating area along cylinder chamber length (hopefully it's a straight line)
    # x_positions = np.arange(0,chamber_length, dx) #position along the chamber length (m)
    # area_values = np.linspace(chamber_area, A_star, len(x_positions)) #local areas (m^2)
    # area_ratios = area_values / A_star #area ratio A/A* (no units)

    # #initializing arrays to store Mach number, heat transfer coefficient, and surface temperature values
    # Mach_array = np.zeros_like(area_ratios) #mach number at each axial position
    # h_array = np.zeros_like(area_ratios) #heat transfer coefficient at each axial position
    # Temp_surface_array = np.zeros_like(area_ratios) #surface temperature at each axial position

    # #The values above were only for the cylindrical chamber section, now adding the nozzle
    # converging_length = 4.464 * c.IN2M #converging section length (m)
    # diverging_length = 1.768 * c.IN2M #diverging section length (m)
    # dx = 0.001 #increments of 1mm

    # #adding new axial positions for converging and diverging sections
    # x_converging = np.arange(-converging_length, 0, dx)
    # x_diverging = np.arange(chamber_length, chamber_length + diverging_length, dx)

    # #adding new area ratio profiles for nozzle
    # #converging section area ratios (from end of cylinder chamber to A*)
    # area_converging = np.linspace(chamber_area, A_star, len(x_converging))
    # area_ratio_converging = area_converging / A_star

    # #diverging section area ratios (from A* to exit area)
    # expansion_ratio = 2.88
    # A_exit = A_star * expansion_ratio
    # area_diverging = np.linspace(A_star, A_exit, len(x_diverging))
    # area_ratio_diverging = area_diverging / A_star

    # #full geometry of the chamber
    # x_total = np.concatenate([x_converging, x_positions, x_diverging])
    # area_ratio_total = np.concatenate([area_ratio_converging, area_ratios, area_ratio_diverging])

    #initializing new arrays to store Mach number, heat transfer coefficient, and surface temperature values for full chamber + nozzle
    Mach_total = np.zeros_like(station_area_ratios) #mach number at each axial position
    h_total = np.zeros_like(station_area_ratios) #heat transfer coefficient at each axial position
    Temp_surface_total = np.zeros_like(station_area_ratios) #surface temperature at each axial position

    initial_guess = 2

    #now calculating Mach number, heat transfer coefficient, and surface temperature at each position along the chamber length
    for station_index, A_ratio in enumerate(station_area_ratios):
        
        M_local = calculating_MachNumber(gamma = cea_results["gamma"], area_ratio_value = A_ratio, initial_guess = initial_guess)
        Mach_total[station_index] = M_local

        #updating initial guess for next iteration
        initial_guess = M_local

        h_local = heat_transfer_coefficient(
            Dt = chamber_diameter, #diameter of chamber (m)
            Rt = D_star,     #radius of throat curve (m)
            Pr = cea_results["c_pran"], #Prandtl number of the combustion gas (n/a)
            gamma = cea_results["gamma"], #specific heat ratio of the combustion gas (n/a)
            c_star = cea_results["c_star"], #characteristic exhaust velocity (m/s)
            T0 = cea_results["c_t"], #stagnation temperature of the combustion gas ((K))
            Cp = cea_results["c_cp"], #specific heat at constant pressure of the combustion gas (J/kg/K)
            P0 = cea_results["c_p"], #chamber pressure (Pascals)
            mu = cea_results["c_visc"], #dynamic viscosity of the combustion gas (Pascal - seconds)
            M = Mach_total[station_index], #Mach number at the local axial point (no units)
            local_Area_ratio = A_ratio #area ratio at the local axial point (no units)

        )
        h_total[station_index] = h_local

        Temp_surface_total[station_index] = temperature_surface_calculation(
            heat_transfer_coefficient_value = h_total[station_index],
            T_infinity = cea_results["c_t"], #chamber temperature (K)
            k = cea_results["c_cond"] #conductivity of the combustion gas in the chamber (W/(m*K))
        )

    
    #printing results

    #printing axial positions vs surface temp plot
    plt.figure()
    plt.plot(station_depths * c.M2IN, Temp_surface_total)
    plt.xlabel("Axial Position Relative to Throat (in) ")
    plt.ylabel("Surface Temperature (K) ")
    plt.title("help")
    plt.grid(True)
    plt.show()

    #printing axial position vs heat transfer coefficient
    plt.figure()
    plt.plot(station_depths * c.M2IN, h_total)
    plt.xlabel("Axial Position Relative to Throat (in) ")
    plt.ylabel("Heat Transfer Coefficient (W/m^2 K) ")
    plt.title("help2")
    plt.grid(True)
    plt.show()

    '''     
    for i in range(len(x_positions)):
        print(f"Axial Position: {x_positions[i]:.3f} m, Mach Number: {Mach_array[i]:.4f}, Heat Transfer Coefficient: {h_array[i]:.2f} W/m^2K, Surface Temperature: {Temp_surface_array[i]:.2f} K")
    '''
    

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
    return float(M_solution[0])


def heat_transfer_coefficient(Dt, Rt, Pr, gamma, c_star, T0, Cp, P0, mu, M, local_Area_ratio):
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


if __name__ == "__main__":
    main()
