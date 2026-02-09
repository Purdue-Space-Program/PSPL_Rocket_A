from CoolProp.CoolProp import PropsSI
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *
import vehicle_parameters as vehicle

def calc_K_sharp_edged_orifice(d, d0):
    beta = d0 / d
    jet_velocity_ratio = 1 + 0.622 * (1 - 0.215 * beta**2 - 0.785 * beta**5)
    K = 0.0696 * (1 - beta**5) * jet_velocity_ratio**2 + (1 - beta**2)**2
    return K

##### Find total area of the orifices needed to sustain a maximum of 30% film #####
max_film_percent = 30 # [%]
m_dot_ipa = vehicle.parameters.core_fuel_mass_flow_rate # [kg/s]
manifold_pressure = vehicle.parameters.chamber_pressure / 0.8 # [Pa] 20% Stiffness to prevent backflow
chamber_pressure = vehicle.parameters.chamber_pressure # [Pa]
dp = manifold_pressure - chamber_pressure # [Pa]
cd = 0.6
max_m_dot = m_dot_ipa * max_film_percent / 100
density = DENSITY_IPA
area_total = max_m_dot / (cd * np.sqrt(2 * density * dp))
print("STEP ONE: Calculating Total Area of Orifices")
print(f"MAX FILM PERCENT: {max_film_percent:.1f}%")
print(f"TOTAL FILM INJECTION AREA: {area_total * M22IN2:.3f} in^2\n")

"""
##### Graph dp vs mdot through adjustable tap-off orifice and film mass flow rate #####

diameter_history = [] # [in]
results = {}
diameter = 0.18 # [in] Smallest Diameter, probably inner diameter of tap off tube

# ENSURE CORRECT DIAMETER UNITS!!!!!!!!!!
while diameter < 0.4:
    area = np.pi * diameter**2 / 4 * IN22M2
    film_percent_history = []
    pressure_drop_history = []
    m_dot_history = []
    p2_history = []
    film_percent = 0
    while film_percent <= 30:
        m_dot = m_dot_ipa * film_percent / 100
        dp = m_dot**2 / (2 * area**2 * cd**2 * density)
        p2 = p1 - dp
        m_dot_history.append(m_dot)
        film_percent_history.append(film_percent)
        pressure_drop_history.append(dp)
        p2_history.append(p2)
        film_percent += 1
    results[diameter] = {
            "film": film_percent_history,
            "mdot": m_dot_history,
            "dp": pressure_drop_history,
            "p2": p2_history
        }
    diameter_history.append(diameter)
    diameter += 0.01 """

def calc_K_sharp_edged_orifice(d, d0):
    beta = d0 / d
    jet_velocity_ratio = 1 + 0.622 * (1 - 0.215 * beta**2 - 0.785 * beta**5)
    K = 0.0696 * (1 - beta**5) * jet_velocity_ratio**2 + (1 - beta**2)**2
    return K

d = 1.5 * IN2M # Diameter of manifold hole 
d0 = 0.05 * IN2M # Diameter of orifice plate
di = 0.18 * IN2M # Inner Diameter of tube before metering orifice
upstream_pressure = 350 * PSI2PA # Pressure before metering orifice, to be changed
dp = np.inf # Initializing pressure drop
film_percent = 0 # Initializing
film_percent_history = []
m_dot_film_history = []
dp_history = []
d0_history = []

##### Find d0 for all mass flow rates for which desired manifold pressure is achieved #####
area_i = np.pi * di**2 / 4 # Inner area of tube before metering orifice

while film_percent <= max_film_percent:
    m_dot_film = m_dot_ipa * film_percent / 100
    line_velocity = m_dot_film / (density * area_i)
    d0 = 0.001 * IN2M
    dp = np.inf
    target_dp = upstream_pressure - manifold_pressure

    while dp > target_dp and d0 < d:
        K0 = calc_K_sharp_edged_orifice(d, d0)
        dp = density * (line_velocity**2 / 2) * K0 * (d**4 / d0**4)
        if dp > target_dp:
            d0 += 0.001 * IN2M
    film_percent_history.append(film_percent)
    m_dot_film_history.append(m_dot_film)
    dp_history.append(dp)
    d0_history.append(d0)
    film_percent += 0.1

plt.plot(np.array(d0_history) * M2IN, film_percent_history)
plt.xlabel("Diameter of Orifice [in]")
plt.ylabel("Film %")
plt.title(f"Diameter of Orifice vs Film %\n[Manifold Diameter = {d * M2IN:.3f}in, Upstream Pressure = {upstream_pressure * PA2PSI:.3f}psi, dp = {dp * PA2PSI:.3f}psi]")
plt.grid()
plt.show()

plt.plot(film_percent_history, np.array(dp_history) * PA2PSI)
plt.xlabel("Film %")
plt.ylabel("Pressure Drop [psi]")
plt.title("Film % vs Pressure Drop")
plt.show()