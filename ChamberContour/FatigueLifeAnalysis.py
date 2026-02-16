import sys
import os
import pandas as pd
sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'utils'))

import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from tqdm import tqdm
from scipy.optimize import root
from scipy.optimize import fsolve
from scipy.optimize import brentq
from scipy.optimize import newton
#from utils.MATERIAL_PROPERTIES import get_material_properties
import Chambercontour as chamber
import vehicle_parameters as vehicle

xchamb, ychamb = chamber.nozzle_contour(vehicle.parameters.chamber_throat_diameter, chamber.expansion_ratio, chamber.Lstar, vehicle.parameters.contraction_ratio, chamber.con_angle, vehicle.parameters.chamber_inner_diameter, chamber.filename)

P_engine_arr = np.array([5e6, 4e6])
P_bulk_arr = np.array([4e6, 3e6])
T_wg_arr = np.array([800, 700])
T_wl_arr = np.array([600, 500])
T_bulk_arr = np.array([400, 300])
material = "carbon steel"
wall_arr = np.array([0.005, 0.005])
width_arr = np.array([0.02, 0.02])
height_arr = np.array([0.02, 0.02])
q_flux_arr = np.array([1e6, 5e5])
T_outer_amb = np.array([300, 300])
T_outer_cool = np.array([280, 280])


def fatigue_eq(n, e, sigma_u, E, strain):

     return 3.5*(sigma_u/E) * n**(-0.12) + (e/n)**(0.6) - strain


def calculate_strains_stresses(P_engine_arr, P_bulk_arr, T_wg_arr, T_wl_arr, T_bulk_arr, material, wall_arr, width_arr, height_arr, q_flux_arr, xchamb, T_outer_amb, T_outer_cool):
    
    print("Running structural analysis...")
    # Get material properties

   
    #axial_stations = len(xchamb)
    axial_stations = len(P_engine_arr) #TEMPPPP

    N_arr = np.array([])
    hotfires_arr = np.array([])
    n_arr = np.array([])
    fires_arr = np.array([])

    for i in tqdm(range(0, axial_stations)):

          # Intialize values
          P_engine = P_engine_arr[i]
          P_bulk = P_bulk_arr[i]
          T_wg = T_wg_arr[i]
          T_wl = T_wl_arr[i]
          T_bulk = T_bulk_arr[i]
          width = width_arr[i]
          height = height_arr[i]
          q_flux = q_flux_arr[i]
          T_outer_a = T_outer_amb[i]
          T_outer_c = T_outer_cool[i]
          t_w = wall_arr[i]

          Dh = (4 * width * height) / (2 * width + 2 * height)
          r = Dh / 2

          # define material temp as average of gas-side and liquid-side temps
          temp = (T_wg + T_wl) / 2

          # Call material properties
          #properties = get_material_properties(material, temp)
          sigma_y = 4.15e+8 #yield strength in Pa
          sigma_u = 5.4e+8 #ultimate strength in Pa
          E = 2e+11 #Young's modulus in GPa
          k = 51.9 #thermal conductivity in W/mK
          a = 0.000012 #thermal expansion coefficient in 1/K
          v = 0.29 #Poisson's ratio
          e = .1 #fracture at elongation

          # Calculate Temperature Gradients
          delta_T1 = T_wg - T_wl
          delta_T2 = ((T_wg + T_wl) / 2) - T_bulk
          delta_T3 = (T_outer_a + T_outer_c) / 2 - T_bulk

          # Calculate thermal stress in unconstrained thin walled tube
          sigma_th1 = (E * a * delta_T1) / (2 * (1 - v))
          # Calculate thermal stress due to initial temperature rise
          sigma_th2 = E * a * delta_T2
          #print(f"Axial Station: {i}")

          sigma_th_axial = sigma_th1 + sigma_th2
         

          # Calculate axial stress due to pressure
          sigma_a_press = (P_bulk * r) / (2 * t_w)

          # Total stresses in axial direction
          sigma_a = sigma_th_axial + sigma_a_press

          # Total strain in axial direction
          strain_tot_a = sigma_a / E

          # Calculate hoop stress
          sigma_h = ((P_bulk - P_engine) * r) / t_w
          sigma_h2 = ((P_bulk - P_engine) / 2) * (width / t_w)**2

          # Total stresses in tangential direction
          sigma_t = sigma_th1 + sigma_h
          # Total strain in tangential direction
          strain_tot_t = sigma_t / E
          # Calculate plastic strain in axial direction
          strain_a_p = strain_tot_a - sigma_y / E
          # Calculate plastic strain in tangential direction
          strain_t_p = strain_tot_t - sigma_y / E
          # Calculate max hot wall temp as to not yield for given heat flux

          # Von mises criterion for plane stress using principal stresses
          J2 = 1/3 * (sigma_t ** 2 - sigma_t * sigma_a + sigma_a ** 2)
          sigma_eff = np.sqrt(3 * J2)
          if sigma_eff > sigma_y:
               # Calculate approximate effective plastic strain
               strain_eff_p = 2 / np.sqrt(3) * np.sqrt(strain_t_p ** 2 + strain_a_p * strain_t_p + strain_a_p ** 2)

          else:
               strain_eff_p = 0

          # Calculate low-cycle Fatigue using eq from SP-8087
          strain_cyclic = strain_eff_p * 2
          N_inverse = (2 * (strain_cyclic - (2 * (0.002 * E)) / E) / e)**2
          N = 1 / N_inverse
          hotfires = N / 4

          # Calculate life using coffin manson relation
          #n = fsolve(fatigue_eq, 70, args=(e, sigma_u, E, strain_eff_p))
          #E = sigma_y / 0.002
          try: 
               n = brentq(fatigue_eq, 1, 1000, args=(e, sigma_u, E, strain_cyclic))
               #n = newton(fatigue_eq, x0=100, args=(e, sigma_y, E, strain_cyclic))
          except:
               n = fsolve(fatigue_eq, 70, args=(e, sigma_u, E, strain_cyclic))
          #print(f"n: {n}")
          fires = n / 4

          N_arr = np.append(N_arr, N)
          hotfires_arr = np.append(hotfires_arr, hotfires)
          n_arr = np.append(n_arr, n)
          fires_arr = np.append(fires_arr, fires)

    return N_arr, hotfires_arr, n_arr, fires_arr



N_arr, hotfires_arr, n_arr, fires_arr = calculate_strains_stresses(P_engine_arr, P_bulk_arr, T_wg_arr, T_wl_arr, T_bulk_arr, material, wall_arr, width_arr, height_arr, q_flux_arr,xchamb, T_outer_amb, T_outer_cool)

print("N_arr:", N_arr)
print("hotfires_arr:", hotfires_arr)
print("n_arr:", n_arr)
print("fires_arr:", fires_arr)