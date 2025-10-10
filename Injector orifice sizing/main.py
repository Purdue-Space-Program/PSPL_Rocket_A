from CoolProp.CoolProp import PropsSI
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *

standard_bits_inch = {
    # Fractional drill sizes (partial list)
    "1/64": 1/64, "1/32": 1/32, "3/64": 3/64, "1/16": 1/16, "5/64": 5/64,
    "3/32": 3/32, "7/64": 7/64, "1/8": 1/8, "9/64": 9/64, "5/32": 5/32,
    "11/64": 11/64, "3/16": 3/16, "13/64": 13/64, "7/32": 7/32, "15/64": 15/64,

    # Number drill sizes (#80–#1, partial list)
    "#80": 0.0135, "#70": 0.0280, "#60": 0.0400, "#50": 0.0700,
    "#40": 0.0980, "#30": 0.1285, "#20": 0.1610, "#10": 0.1935, "#1": 0.2280,
}

def closest_bit_size(target_diameter):
    lowest_absolute_error_bit_size = None

    for key, value in standard_bits_inch.items():
        absolute_error_bit_size = abs(value - target_diameter)
        if (lowest_absolute_error_bit_size is None) or (absolute_error_bit_size < lowest_absolute_error_bit_size):
            lowest_absolute_error_bit_size = absolute_error_bit_size
            closest_bit_name = key
            closest_bit_diameter = value

    return closest_bit_diameter * IN2M, lowest_absolute_error_bit_size * IN2M, closest_bit_name

# Inputs
C_D_lox = 0.6 # Discharge coefficient of oxidizer orifice 
C_D_ipa = 0.6 # Discharge coefficient of ipa orifice
N_lox_max = 100 # Max amount of orifices allowed, to be changed
N_lox_min = 10 # Min amount of orifices allowed, to be changed
of_ratio = 1 # from Vehicle Parameters page
D_c = 4.5 * IN2M # Diameter of chamber (might be changed)
m_dot = 3.56 * LB2KG # [Kg/s]
temp_lox = 90 # [K]
pressure_lox = 250 * PSI2PA # [Pa]
temp_ipa = 290 # [K]
pressure_ipa = 250 * PSI2PA # [Pa]
rho_lox = PropsSI('D', 'T', temp_lox, 'P', pressure_lox, 'Oxygen')
rho_ipa = DENSITY_IPA # [kg/m^3] CoolProp doesn't have IPA, assuming constant
pressure_chamber = 150 * PSI2PA # [Pa]
pressure_upstream_injector = pressure_chamber / 0.8 # [Pa]
pressure_drop = pressure_upstream_injector - pressure_chamber # [Pa]
m_dot_ipa = m_dot / (1 + of_ratio) # [Kg/s]
m_dot_ideal_lox = m_dot - m_dot_ipa # [Kg/s]
D_s = D_c / 5 # Diameter of pintle shaft
R_s = D_s / 2 # Radius of pintle shaft
skip_dist = D_s # This means that the skip distance ratio = 1. This is a good rule of thumb.

total_target_area_orifice_lox = m_dot_ideal_lox / (C_D_lox * np.sqrt(2 * rho_lox * pressure_drop))

total_area_orifice_ipa = m_dot_ipa / (C_D_ipa * np.sqrt(2 * rho_ipa * pressure_drop))
velocity_ipa = m_dot_ipa / (rho_ipa * total_area_orifice_ipa)
annular_thickness = np.sqrt((total_area_orifice_ipa + np.pi * R_s**2) / np.pi) - R_s

# Constants and initializing variables
g = 9.81 # [m/s^2]
N_rows = 1
tmr_ideal = 1.3 # Optimal total momentum ratio, as stated in Fundamental Combustion Characteristics of Ethanol/Liquid Oxygen Rocket Engine Combustor with Planar Pintle-type Injector"
tmr_allowable_percent_error = 0.4 # percent error allowed in pintle's TMR
allowable_percent_error_m_dot_lox = 0.006 # percent error allowed in LOx mass flow rate
good_enough_found = 0

N_lox_array = []
value_1_array = []
value_2_array = []

for N_lox in range(N_lox_min, N_lox_max, 2):
    # print(f"N_lox: {N_lox}")
    
    # consider that fact that bits in real life are not the exact diameter we want
    D_ideal_lox_orifice = 2 * np.sqrt(total_target_area_orifice_lox / (np.pi * N_lox))
    # print(f"D_target_lox_orifice: {D_target_lox_orifice * M2IN:.2f}")
    closest_bit_diameter, absolute_error_bit_size, closest_bit_name = closest_bit_size(D_ideal_lox_orifice * M2IN)

    D_real_lox_orifice = closest_bit_diameter
    percent_error_bit_size = abs(absolute_error_bit_size) / D_ideal_lox_orifice
    
    total_area_real_orifice_lox = ((D_real_lox_orifice/2)**2) * (np.pi * N_lox)
    # print(f"total_area_real_orifice_lox: {total_area_real_orifice_lox * M2IN:.2}")
    m_dot_real_lox = total_area_real_orifice_lox * (C_D_lox * np.sqrt(2 * rho_lox * pressure_drop))
    
    percent_error_m_dot_lox = abs(m_dot_real_lox - m_dot_ideal_lox) / m_dot_ideal_lox
    
    # print(f"area_orifice_lox: {total_area_real_orifice_lox:.2f}")
    # print(f"area_orifice_ipa: {total_area_orifice_ipa:.2f}")
    
    velocity_lox = m_dot_real_lox / (rho_lox * total_area_real_orifice_lox)
    
    
    # print(f"velocity_lox: {velocity_lox:.2f}")
    # print(f"velocity_ipa: {velocity_ipa:.2f}")
    tmr_real = (m_dot_real_lox * velocity_lox) / (m_dot_ipa * velocity_ipa)
    # print(f"tmr_real: {tmr_real:.5f}")

    
    bf = N_lox * D_real_lox_orifice / (np.pi * D_s)
    N_rows = (N_lox * D_real_lox_orifice // (np.pi * D_s)) + 1
    
    lmr = tmr_real/bf
    half_angle = 0.7 * np.degrees(np.arctan(2 * lmr)) # from: 
                                                                    # Blakely, J., Freeberg, J., and Hogge, J., “Spray Cone Formation from
                                                                    # Pintle-Type Injector Systems in Liquid Rocket Engines,” AIAA SciTech
                                                                    # 2019 Forum, AIAA Paper 2019-0152, 2019.
                                                                    # https://doi.org/10.2514/6.2019-0152
    tmr_percent_error = abs(tmr_real - tmr_ideal) / tmr_ideal
    # print(f"tmr_error: {tmr_error:.5f}")
    # print(tmr_error < tmr_allowable_error)
    
    # print(f"tmr_error < tmr_allowable_error: {tmr_error <= tmr_allowable_error}")
    
    
    N_lox_array.append(N_lox)
    value_1_array.append(tmr_percent_error)
    value_2_array.append(percent_error_m_dot_lox)
    
    # only update orifice sizes if a good enough one has not been found
    if (tmr_percent_error <= tmr_allowable_percent_error) and (percent_error_m_dot_lox <= allowable_percent_error_m_dot_lox) and (good_enough_found == 0):
        good_enough_found = 1
        best_values = {
            "Number of holes": N_lox,
            "Number of rows": N_rows,
            "Closest Bit Name": closest_bit_name,
            
            "Ideal diameter of LOx orifice holes [in]": D_ideal_lox_orifice * M2IN,
            "Real diameter of LOx orifice holes [in]": D_real_lox_orifice * M2IN,
            "Bit size error [%]": percent_error_bit_size * 100,
            
            "Ideal LOx mass flow rate [kg/s]" : m_dot_ideal_lox,
            "Real LOx mass flow rate [kg/s]" : m_dot_real_lox,
            "LOx mass flow rate error [%]": percent_error_m_dot_lox * 100,
            
            "Blockage Factor": bf,
            "TMR": tmr_real,
            "Optimal TMR": tmr_ideal,
            "TMR Error from optimal [%]": tmr_percent_error * 100,
            "LMR": lmr,
            "half_angle": half_angle,
            "velocity_lox [m/s]": velocity_lox,
            "velocity_ipa [m/s]": velocity_ipa,
            "area_orifice_ipa [in^2]": total_area_orifice_ipa * M22IN2,
            "area_orifice_lox [in^2]": total_area_real_orifice_lox * M22IN2,
            }
    

if good_enough_found == 1:
    for key, value in best_values.items():
        if isinstance(value, str):
            print(f"{key}: {value}")
        else:
            print(f"{key}: {value:.3f}")
    # show where chosen orifice size is
    plt.axvline(best_values["Number of holes"], color='g', linestyle='--', label='Chosen number of holes')
else:
    print("NO GOOD PINTLE FOUND (ERRORS TOO HIGH)")


plt.plot(N_lox_array, value_1_array, color='b', label="tmr_percent_error")
plt.plot(N_lox_array, value_2_array, color='black', label="percent_error_m_dot_lox")


# show where allowable errors are
plt.axhline(tmr_allowable_percent_error, color='orange', linestyle='--', label='Max allowable percent error in TMR')
plt.axhline(allowable_percent_error_m_dot_lox, color='r', linestyle=':', label='Max allowable percent error in LOx mass flow rate')

plt.grid(True, which='both', axis='x', color='lightgray', linestyle='-', linewidth=0.5)
plt.xticks(np.arange(min(N_lox_array), max(N_lox_array) + 1, 2))
plt.legend()
plt.xlabel("Number of LOx holes")
# plt.ylabel("LMR")
plt.show()
