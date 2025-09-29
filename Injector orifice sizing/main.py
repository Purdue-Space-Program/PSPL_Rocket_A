from CoolProp.CoolProp import PropsSI
import numpy as np
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *

standard_bits_inch = {
    # Fractional drill sizes (partial list)
    "1/64": 1/64, "1/32": 1/32, "3/64": 3/64, "1/16": 1/16, "5/64": 5/64,
    "3/32": 3/32, "7/64": 7/32, "1/8": 1/8, "9/64": 9/64, "5/32": 5/32,
    "11/64": 11/64, "3/16": 3/16, "13/64": 13/64, "7/32": 7/32, "15/64": 15/64,
    "1/4": 1/4, "5/16": 5/16, "3/8": 3/8, "7/16": 7/16, "1/2": 1/2,

    # Number drill sizes (#80–#1, partial list)
    "#80": 0.0135, "#70": 0.0280, "#60": 0.0400, "#50": 0.0700,
    "#40": 0.0980, "#30": 0.1285, "#20": 0.1610, "#10": 0.1935, "#1": 0.2280,

    # Letter drill sizes (A–Z, partial list)
    "A": 0.234, "B": 0.238, "C": 0.242, "D": 0.246, "E": 0.250,
    "F": 0.257, "G": 0.261, "H": 0.266, "I": 0.272, "J": 0.277,
    "K": 0.281, "L": 0.290, "M": 0.295, "N": 0.302, "O": 0.316,
    "P": 0.323, "Q": 0.332, "R": 0.339, "S": 0.348, "T": 0.358,
    "U": 0.368, "V": 0.377, "W": 0.386, "X": 0.397, "Y": 0.404, "Z": 0.413,
}
g = 9.81 # [m/s^2]
C_D_lox = 0.795 # Discharge coefficient of oxidizer orifice 
C_D_ipa = 0.948 # Discharge coefficient of ipa orifice
m_dot = 3.56 * LB2KG # [Kg/s]
temp_lox = 90 # [K]
pressure_lox = 250 * PSI2PA # [Pa]
temp_ipa = 290 # [K]
pressure_ipa = 250 * PSI2PA # [Pa]
rho_lox = PropsSI('D', 'T', temp_lox, 'P', pressure_lox, 'Oxygen')
rho_ipa = DENSITY_IPA # [kg/m^3] CoolProp doesn't have IPA, assuming constant
pressure_chamber = 150 * PSI2PA # [Pa]
pressure_upstream_injector = pressure_chamber / 0.8 # [Pa]
of_ratio = 1 # from Vehicle Parameters page
m_dot_ipa = m_dot / (1 + of_ratio) # [Kg/s]
m_dot_lox = m_dot - m_dot_ipa # [Kg/s]
D_c = 4.5 * IN2M # Diameter of chamber (might be changed)
D_s = D_c / 5 # Diameter of pintle shaft
R_s = D_s / 2 # Radius of pintle shaft
skip_dist = D_s # This means that the skip distance ratio = 1. This is a good rule of thumb.
pressure_drop = pressure_upstream_injector - pressure_chamber # [Pa]
total_area_orifice_lox = m_dot_lox / (C_D_lox * np.sqrt(2 * rho_lox * pressure_drop))
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

def closest_bit_size(d):
    closest_bit = None
    min_distance = None
    for i in standard_bits_inch.items():
        distance = abs(i[1] - d)
        if min_distance is None or distance < min_distance:
            min_distance = distance
            closest_bit = i
    return d, closest_bit, min_distance

while N_lox <= N_lox_max:
    D_lox_orifice = 2 * np.sqrt(total_area_orifice_lox / (np.pi * N_lox))
    D_lox_orifice_inches, closest_bit, min_distance = closest_bit_size(D_lox_orifice * M2IN)
    print(f"Diameter: {D_lox_orifice:.5f} m = {D_lox_orifice_inches:.5f} in")
    print(f"Closest drill bit: {closest_bit}")
    print(f"Difference: {min_distance:.5f} in")
    D_lox_orifice_real = closest_bit[1] * IN2M
    area_orifice_lox = total_area_orifice_lox / N_lox
    annular_thickness = np.pi * rho_lox * D_lox_orifice_real / (4 * rho_ipa * of_ratio**2)
    area_orifice_ipa = np.pi * (R_s + annular_thickness)**2 - np.pi * R_s**2
    velocity_ipa = m_dot_ipa / (rho_ipa * area_orifice_ipa)
    velocity_lox = m_dot_lox / (rho_lox * area_orifice_lox)
    tmr = (m_dot_lox * velocity_lox) / (m_dot_ipa * velocity_ipa)
    bf = N_lox * D_lox_orifice_real / (np.pi * D_s)
    if bf > 1:
        N_rows += 1
    lmr = tmr/bf
    half_angle = 0.7 * np.arctan(2 * lmr)
    print(f"TMR: {tmr}")
    

    N_lox += 1
    break

print(f"LOx diameter: {D_lox_orifice:.3} meters")
print(f"LOx diameter: {D_lox_orifice * M2IN:.3} inches")
