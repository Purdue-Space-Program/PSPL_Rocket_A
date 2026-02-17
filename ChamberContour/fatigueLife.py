import numpy as np
from scipy.optimize import brentq
from tqdm import tqdm
import Chambercontour as chamber
import matplotlib.pyplot as plt

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
import vehicle_parameters as vehicle
import pandas as pd



def pressure_from_mach(P0, Mach_arr, gamma=1.22):
    return P0 * (1 + (gamma - 1)/2 * Mach_arr**2)**(-gamma/(gamma - 1))

def calculating_MachNumber(gamma, area_ratio_value, initial_guess, branch):

    def f(M):
        Mach_function_part1 = ((gamma + 1)/2)**(-(gamma + 1)/(2*(gamma-1)))
        Mach_function_part2 = (1/M) * ((1 + (((gamma - 1)/2) * M**2)) ** ((gamma + 1)/(2*(gamma-1))))
        Mach_function = (Mach_function_part1 * Mach_function_part2) - area_ratio_value
        return Mach_function

    guess = initial_guess

    
    if branch == "subsonic":
        guess = min(0.8, max(0.05, guess))

    else:
        guess = max(1.1, guess)
    
    M_solution = fsolve(f, guess)
    return float(M_solution[0])


def cycles_to_failure_coffin_manson(strain_eff_p, eps_f, c=-0.6): #for low cycle fatigue
    delta_eps_p = 2.0 * strain_eff_p
    return 0.5 * (delta_eps_p / (2 * eps_f))**(1.0 / c)


def calculate_strains_stresses_contour(xchamb, ychamb, P_engine_arr, T_wg_arr, wall, T_outer_amb):
    
    axial_stations = len(xchamb)

    sigma_y = 68.94757 # yield strength (MPa)
    E = 200000 # Young's modulus(MPa)
    alpha = 12e-6 # thermal expansion coeff
    nu = 0.3 # Poisson's ratio
    eps_f = 0.50 # fatigue ductility coeff
    c = -0.6

    Nf_arr = np.zeros(axial_stations)
    strain_arr =  np.zeros(axial_stations)

    for i in range(axial_stations):

        r_inner = ychamb[i]
        t_w = wall
        r_mean = r_inner + t_w/2

        # Thermal stress
        delta_T = T_wg_arr[i] - T_outer_amb[i]
        sigma_th = (E * alpha * delta_T) / (1 - nu)

        # Pressure stresses
        sigma_hoop = P_engine_arr[i] * r_mean / t_w
        sigma_axial = P_engine_arr[i] * r_mean / (2 * t_w)

        # Total stresses
        sigma_theta = sigma_th + sigma_hoop
        sigma_x = sigma_th + sigma_axial

        # Elastic strains
        eps_theta = sigma_theta / E #hooke's law
        eps_x = sigma_x / E

        # Plastic strains
        eps_theta_p = max(eps_theta - sigma_y/E, 0)
        eps_x_p = max(eps_x - sigma_y/E, 0)

        eps_eff_p = (2/np.sqrt(3)* np.sqrt(eps_theta_p**2 + eps_x_p**2 + eps_theta_p * eps_x_p)) #combines theta and x strains
        
        strain_arr[i] = eps_eff_p


        Nf_arr[i] = cycles_to_failure_coffin_manson(eps_eff_p, eps_f, c)

    plt.figure(1)
    plt.plot(Nf_arr, strain_arr)
    plt.xlabel("cycles")
    plt.ylabel("strain")
    plt.title("strain v cycles")
    plt.show()


    return Nf_arr



xchamb, ychamb = chamber.nozzle_contour(vehicle.parameters.chamber_throat_diameter, chamber.expansion_ratio, chamber.Lstar, vehicle.parameters.contraction_ratio, chamber.con_angle, vehicle.parameters.chamber_inner_diameter, chamber.filename)
xchamb = c.IN2M * xchamb
ychamb = c.IN2M * ychamb

# Area from contour
A_arr = np.pi * ychamb**2
throat_index = np.argmin(A_arr)
At = A_arr[throat_index]

Mach_arr = np.zeros_like(A_arr)

for i, A in enumerate(A_arr):
    Ar = A / At

    if np.isclose(Ar, 1.0, atol=1e-6):
        Mach_arr[i] = 1.0

    elif i < throat_index:
        Mach_arr[i] = calculating_MachNumber(max(Ar, 1.0), supersonic=False)

    else:
        Mach_arr[i] = calculating_MachNumber(max(Ar, 1.0), supersonic=True)


Pc = 1.724  # chamber pressure (MPa)
gamma = 1.22

P_engine_arr = pressure_from_mach(Pc, Mach_arr, gamma)
np.savetxt("chamber_pressure.csv", P_engine_arr, delimiter=',')

T_wg_arr = pd.read_csv('ChamberContour/chamber_temp.csv').to_numpy().ravel() # wall gas side temp [K]
wall = 0.00635  # wall thickness [m]
T_outer_amb = np.linspace(300, 300, len(xchamb))  # outer ambient temp [K]

Nf_arr = calculate_strains_stresses_contour(xchamb,ychamb,P_engine_arr,T_wg_arr,wall,T_outer_amb)

print("\nMinimum chamber life (cycles):", np.min(Nf_arr))
print("Critical location x (m):", xchamb[np.argmin(Nf_arr)])