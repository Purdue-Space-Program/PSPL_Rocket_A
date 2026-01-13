from CoolProp.CoolProp import PropsSI
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *
import vehicle_parameters as vehicle

standard_bits_inch = {
    # Fractional drill sizes (partial list)
    "1/64": 1/64, "1/32": 1/32, "3/64": 3/64, # may not have these
    "1/16": 1/16, "5/64": 5/64, "3/32": 3/32, # for sure have these
    "7/64": 7/64, "1/8": 1/8, # idk, probably too big

    # Number drill sizes (#80–#1, partial list)
    "#80": 0.0135, "#70": 0.0280, "#60": 0.0400, "#59": 0.0410, "#58": 0.0420,
    "#57": 0.0430, "#56": 0.0465, "#55": 0.0520, "#54": 0.0550, "#53": 0.0595,
    "#52": 0.0635, "#51": 0.0670, "#50": 0.0700, "#49": 0.0730, "#48": 0.0760,
    "#47": 0.0785, "#46": 0.0810, "#45": 0.0820, "#44": 0.0860, "#43": 0.0890,
    "#42": 0.0935, "#41": 0.0960, "#40": 0.0980, "#39": 0.0995, "#38": 0.1015,
}

standard_bits_metric = {
    # Metric bit sizes (partial list available at BIDC)
    "1 mm":   0.0010, "1.1 mm": 0.0011, "1.2 mm": 0.0012, "1.3 mm": 0.0013,
    "1.4 mm": 0.0014, "1.5 mm": 0.0015, "1.6 mm": 0.0016, "1.7 mm": 0.0017,
    "1.8 mm": 0.0018, "1.9 mm": 0.0019, "2 mm":   0.0020, "2.1 mm": 0.0021,
    "2.2 mm": 0.0022, "2.3 mm": 0.0023, "2.4 mm": 0.0024, "2.5 mm": 0.0025,
    "film exception": 10000000000000000000000000000000000000000000,
}

remove_list = ["1/64", "1/32", "3/64", "#80", "#70"]

for key, value in standard_bits_inch.items():
    if key not in remove_list:
        standard_bits_metric[key] = value * IN2M 


def closest_bit_size(target_diameter, already_used_bit_size):
    lowest_absolute_error_bit_size = None
    
    fuck = standard_bits_metric.copy()
    fuck.pop(already_used_bit_size)

    for key, value in fuck.items():
        absolute_error_bit_size = abs(value - target_diameter)
        if (lowest_absolute_error_bit_size is None) or (absolute_error_bit_size < lowest_absolute_error_bit_size):
            lowest_absolute_error_bit_size = absolute_error_bit_size
            closest_bit_name = key
            closest_bit_diameter = value

    return closest_bit_diameter, lowest_absolute_error_bit_size, closest_bit_name

def CalculateAreaFromMassFlowRate(m_dot, C_D, rho, pressure_drop):
    total_area = m_dot / (C_D * np.sqrt(2 * rho * pressure_drop))
    return total_area

def CalculateAreaFromHoles(hole_diameter, number_holes):
    hole_radius = hole_diameter/2
    area_per_hole = np.pi * (hole_radius)**2
    total_area = area_per_hole * number_holes

    return total_area

def CalculateMassFlowRate(area, Cd, rho, pressure_drop):
    m_dot = area * (Cd * np.sqrt(2 * rho * pressure_drop))
    return m_dot
def CalculateIdealHoleDiameter(area, number_holes):
    area_per_hole = area/number_holes
    ideal_hole_radius = np.sqrt(area_per_hole / np.pi)
    ideal_hole_diameter = 2 * ideal_hole_radius
    return ideal_hole_diameter

# Inputs
C_D_lox = 0.6 # Discharge coefficient of oxidizer orifice 
C_D_ipa = 0.6 # Discharge coefficient of ipa orifice
N_bottom_max = 40 # Max amount of orifices in bottom row allowed, to be changed (we will have 2 rows where the bottom row will have holes with smaller diameter)
N_bottom_min = 8 # Min amount of orifices in bottom row allowed, to be changed
of_ratio = vehicle.parameters.OF_ratio # from Vehicle Parameters page
D_c = vehicle.parameters.chamber_inner_diameter # Diameter of chamber (might be changed)
m_dot = vehicle.parameters.total_mass_flow_rate # [Kg/s]
temp_lox = 90 # [K]
pressure_lox = vehicle.parameters.tank_pressure # [Pa]
temp_ipa = 290 # [K]
pressure_ipa = vehicle.parameters.tank_pressure # [Pa]
rho_lox = PropsSI('D', 'T', temp_lox, 'P', pressure_lox, 'Oxygen')
rho_ipa = DENSITY_IPA # [kg/m^3] CoolProp doesn't have IPA, assuming constant
pressure_chamber = vehicle.parameters.chamber_pressure # [Pa]
pressure_upstream_injector = pressure_chamber / 0.8 # [Pa]
desired_pressure_drop = pressure_upstream_injector - pressure_chamber # [Pa]
m_dot_ipa = m_dot / (1 + of_ratio) # [Kg/s]
m_dot_ideal_lox = m_dot - m_dot_ipa # [Kg/s]
D_s = D_c / 5 # Diameter of pintle shaft
R_s = D_s / 2 # Radius of pintle shaft
skip_dist = D_s # This means that the skip distance ratio = 1. This is a good rule of thumb.
total_target_area_orifice_lox = CalculateAreaFromMassFlowRate(m_dot_ideal_lox, C_D_lox, rho_lox, desired_pressure_drop)
total_area_orifice_ipa = CalculateAreaFromMassFlowRate(m_dot_ipa, C_D_ipa, rho_ipa, desired_pressure_drop)
velocity_ipa = m_dot_ipa / (rho_ipa * total_area_orifice_ipa)
annular_thickness = np.sqrt((total_area_orifice_ipa + np.pi * R_s**2) / np.pi) - R_s

# Constants and initializing variables
g = 9.81 # [m/s^2]
N_rows = 1
tmr_ideal = 1.3 # Optimal total momentum ratio, as stated in Fundamental Combustion Characteristics of Ethanol/Liquid Oxygen Rocket Engine Combustor with Planar Pintle-type Injector"
tmr_allowable_percent_error = 40 # percent error allowed in pintle's TMR
allowable_percent_error_m_dot_lox = 1 # percent error allowed in LOx mass flow rate
good_enough_found = 0

N_lox_array = []
value_1_array = []
value_2_array = []
minerr = np.inf
D_real_lox_orifice_top = "#55"
#N_top = 10
for N_top in range(20, 30, 2):
    for N_bottom in range(N_top, N_bottom_max + 1, 2):
        N_lox_array.append(N_bottom)
        N_lox = N_top + N_bottom
        total_area_real_top_orifice_lox = CalculateAreaFromHoles(standard_bits_metric[D_real_lox_orifice_top], N_top)
        total_target_area_bottom_orifice_lox = total_target_area_orifice_lox - total_area_real_top_orifice_lox
        D_ideal_lox_orifice_bottom = CalculateIdealHoleDiameter(total_target_area_bottom_orifice_lox, N_bottom) 
        D_real_lox_orifice_bottom, lowest_absolute_error_bit_size_bottom, closest_bit_name_bottom = closest_bit_size(D_ideal_lox_orifice_bottom, D_real_lox_orifice_top)
        percent_error_bit_size_bottom = (abs(lowest_absolute_error_bit_size_bottom) / D_ideal_lox_orifice_bottom) * 100
        total_area_real_bottom_orifice_lox = CalculateAreaFromHoles(D_real_lox_orifice_bottom, N_bottom)
        total_area_real_orifice_lox = total_area_real_bottom_orifice_lox + total_area_real_top_orifice_lox
        err = (abs(total_area_real_orifice_lox - total_target_area_orifice_lox) / total_target_area_orifice_lox) * 100
        m_dot_real_lox = CalculateMassFlowRate(total_area_real_orifice_lox, C_D_lox, rho_lox, desired_pressure_drop)
        percent_error_m_dot_lox = (abs(m_dot_real_lox - m_dot_ideal_lox) / m_dot_ideal_lox) * 100
        velocity_lox = m_dot_real_lox / (rho_lox * total_area_real_orifice_lox)
        tmr_real = (m_dot_real_lox * velocity_lox) / (m_dot_ipa * velocity_ipa)
        bf = (N_top * standard_bits_metric[D_real_lox_orifice_top] + N_bottom * D_real_lox_orifice_bottom) / (np.pi * D_s)
        lmr = tmr_real/bf
        half_angle = 0.7 * np.degrees(np.arctan(2 * lmr)) # from: 
                                                                    # Blakely, J., Freeberg, J., and Hogge, J., “Spray Cone Formation from
                                                                    # Pintle-Type Injector Systems in Liquid Rocket Engines,” AIAA SciTech
                                                                    # 2019 Forum, AIAA Paper 2019-0152, 2019.
                                                                    # https://doi.org/10.2514/6.2019-0152
        tmr_percent_error = (abs(tmr_real - tmr_ideal) / tmr_ideal) * 100
        value_1_array.append(tmr_percent_error)
        value_2_array.append(percent_error_m_dot_lox)
        # only update orifice sizes if a good enough one has not been found
        if (tmr_percent_error <= tmr_allowable_percent_error) and (percent_error_m_dot_lox <= allowable_percent_error_m_dot_lox) and (good_enough_found == 0):
            good_enough_found = 1
            best_values = {
                "Number of holes": N_bottom + N_top,
                "Number of holes in top row": N_top,
                "Number of holes in bottom row": N_bottom,
                "Skip distance [m]": skip_dist,
                "Annular Thickness [in]": annular_thickness * M2IN,
                "Real diameter of top LOx orifice holes [in]": standard_bits_metric[D_real_lox_orifice_top] * M2IN,
                "Bit name (top orifices)": D_real_lox_orifice_top,
                "Ideal diameter of bottom LOx orifice holes [in]": D_ideal_lox_orifice_bottom * M2IN,
                "Real diameter of bottom LOx orifice holes [in]": D_real_lox_orifice_bottom * M2IN,
                "Closest Bit Name (bottom orifices)": closest_bit_name_bottom,
                "Bit size error of bottom orifice holes [%]": percent_error_bit_size_bottom,
                
                "Ideal LOx mass flow rate [kg/s]" : m_dot_ideal_lox,
                "Real LOx mass flow rate [kg/s]" : m_dot_real_lox,
                "LOx mass flow rate error [%]": percent_error_m_dot_lox,
                
                "Blockage Factor": bf,
                "TMR": tmr_real,
                "Optimal TMR": tmr_ideal,
                "TMR Error from optimal [%]": tmr_percent_error,
                "LMR": lmr,
                "half_angle": half_angle,
                "velocity_lox [m/s]": velocity_lox,
                "velocity_ipa [m/s]": velocity_ipa,
                "area_orifice_ipa [in^2]": total_area_orifice_ipa * M22IN2,
                "area_orifice_lox [in^2]": total_area_real_orifice_lox * M22IN2,
                "top area ratio": total_area_real_top_orifice_lox / total_area_real_orifice_lox,
                "error in area [%]": err
            }


if good_enough_found == 1:
    for key, value in best_values.items():
        if isinstance(value, str):
            print(f"{key}: {value}")
        else:
            print(f"{key}: {value:.5f}")
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
plt.xlabel("Number of Bottom LOx holes")
# plt.ylabel("LMR")
plt.show()

########### Film cooling orifices sizing ##########

C_D_film = 0.6
m_dot_ideal_film = 0.25 * m_dot_ipa / 0.75
#print(m_dot_ideal_film)
total_target_area_orifice_film = CalculateAreaFromMassFlowRate(m_dot_ideal_film, C_D_film, rho_ipa, desired_pressure_drop)
N_film_max = 100
N_film_min = 1
allowable_percent_error_m_dot_film = 3 # percent error allowed in Film mass flow rate
x_num = []
y_err = []

for N_film in range(N_film_min, N_film_max + 1):
    D_ideal_orifice_film = CalculateIdealHoleDiameter(total_target_area_orifice_film, N_film)
    closest_bit_diameter, absolute_error_bit_size_film, closest_bit_name = closest_bit_size(D_ideal_orifice_film, "film exception")
    D_real_orifice_film = closest_bit_diameter
    
    total_area_real_orifice_film = CalculateAreaFromHoles(D_real_orifice_film, N_film)
    m_dot_real_film = CalculateMassFlowRate(total_area_real_orifice_film, C_D_film, rho_ipa, desired_pressure_drop)
    
    percent_error_bit_size_film = (abs(absolute_error_bit_size_film) / D_ideal_orifice_film) * 100
    percent_error_m_dot_film = (abs(m_dot_real_film - m_dot_ideal_film) / m_dot_ideal_film) * 100
    
    if percent_error_m_dot_film <= allowable_percent_error_m_dot_film:
        print(f"\n------------ FILM ------------")
        print(f"Number of film orifices: {N_film}")
        print(f"Diameter of target film hole [in]: {D_ideal_orifice_film * M2IN:.6f}")
        print(f"Diameter of film hole [m]: {D_real_orifice_film:.6f}")
        print(f"Diameter of film hole [in]: {D_real_orifice_film * M2IN:.6f}")
        print(f"Bit size(in): {closest_bit_name}")
        print(f"% Error Bit Size: {percent_error_bit_size_film}")
        print(f"Error percent of mass flow rate: {percent_error_m_dot_film}")
    x_num.append(N_film)
    y_err.append(percent_error_m_dot_film)
    #print(percent_error_m_dot_film)
plt.plot(x_num, y_err)
plt.xlabel("Number of holes")
plt.ylabel("Percentage error in mass flow rate")
plt.title("[Film Cooling] Number of holes vs m_dot error comparision")
plt.axhline(allowable_percent_error_m_dot_film, color='g', linestyle='--', label='Chosen number of holes')
plt.show()