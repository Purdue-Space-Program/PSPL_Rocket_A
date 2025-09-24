from CoolProp.CoolProp import PropsSI
import numpy as np
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *

g = 9.81 # [m/s^2]
C_D_lox = 0.795 # Discharge coefficient of oxidizer orifice 
C_D_ipa = 0.948 # Discharge coefficient of ipa orifice
m_dot = 3.05 * LB2KG # [Kg/s]
temp_lox = 90 # [K]
pressure_lox = 250 * PSI2PA # [Pa]
temp_ipa = 290 # [K]
pressure_ipa = 250 * PSI2PA # [Pa]
rho_lox = PropsSI('D', 'T', temp_lox, 'P', pressure_lox, 'Oxygen')
rho_ipa = DENSITY_IPA # [Pa] CoolProp doesn't have IPA, assuming constant
pressure_chamber = 190 * PSI2PA # [Pa]
pressure_upstream_injector = pressure_chamber / 0.8 # [Pa]
of_ratio = 1 # from Vehicle Parameters page
m_dot_ipa = m_dot / (1 + of_ratio) # [Kg/s]
m_dot_lox = m_dot - m_dot_ipa # [Kg/s]
D_c = 4.5 * IN2M # Diameter of chamber (might be changed)
D_s = D_c / 5 # Diameter of pintle shaft
R_s = D_s / 2 # Radius of pintle shaft
skip_dist = D_s # This means that the skip distance ratio = 1. This is a good rule of thumb.
pressure_drop = pressure_upstream_injector - pressure_chamber# [Pa]
total_area_orifice_lox = m_dot_lox / (C_D_lox * np.sqrt(2 * g * rho_lox * pressure_drop))
N_lox_max = 120 # Max amount of orifices allowed, to be changed
N_lox = 1
N_rows = 1
annular_thickness = 0 # [m]
area_orifice_ipa = 0 # [m^2]
area_orifice_lox = 0 # [m^2]
D_lox_orifice = 0 # [m]
velocity_ipa = 0 # [m/s]
velocity_lox = 0 # [m/s]
tmr_real = 0 # Real total momentum ratio
tmr_optimal = 1 # Optimal total momentum ratio 
bf = 0 # Blockage factor
lmr = 0 # [Local momentum ratio]
half_angle = 0 # [radians]

while N_lox <= N_lox_max:
    D_lox_orifice = 2 * np.sqrt(total_area_orifice_lox / (np.pi * N_lox))
    area_orifice_lox = total_area_orifice_lox / N_lox
    annular_thickness = np.pi * rho_lox * D_lox_orifice / (4 * rho_ipa * of_ratio**2)
    area_orifice_ipa = np.pi * (R_s + annular_thickness)**2 - np.pi * R_s**2
    velocity_ipa = m_dot_ipa / (rho_ipa * area_orifice_ipa)
    velocity_lox = m_dot_lox / (rho_lox * area_orifice_lox)
    tmr = (m_dot_lox * velocity_lox) / (m_dot_ipa * velocity_ipa)
    bf = N_lox * D_lox_orifice / (np.pi * D_s)
    if bf > 1:
        N_rows += 1
    lmr = tmr/bf
    half_angle = 0.7 * np.arctan(2 * lmr)
    print(tmr)

    N_lox += 1
    break

print(D_lox_orifice)