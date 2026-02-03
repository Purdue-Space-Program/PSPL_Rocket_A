from CoolProp.CoolProp import PropsSI
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *
import vehicle_parameters as vehicle

##### Find total area of the orifices needed to sustain a maximum of 30% film #####
max_film_percent = 30 # [%]
m_dot_ipa = vehicle.parameters.fuel_mass_flow_rate # [kg/s]
p1 = 310 * PSI2PA # [Pa]
p2 = vehicle.parameters.chamber_pressure # [Pa]
dp = p1 - p2 # [Pa]
cd = 0.6
max_m_dot = m_dot_ipa * max_film_percent / 100
density = DENSITY_IPA
area_total = max_m_dot / (cd * np.sqrt(2 * density * dp))
print("STEP ONE: Calculating Total Area of Orifices")
print(f"MAX FILM PERCENT: {max_film_percent:.1f}%")
print(f"TOTAL FILM INJECTION AREA: {area_total * M22IN2:.3f} in^2\n")

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
    diameter += 0.01

print("STEP 2: Graphing Pressure Drop Through Orifice with Film Mass Flow Rate")

for d in diameter_history:
    plt.plot(
        results[d]["film"],
        np.array(results[d]["dp"]) * PA2PSI,
        label=f"d = {d:.2f} in"
    )
plt.xlabel("Film Percentage (%)")
plt.ylabel("Pressure Drop (PSI)")
plt.title("Pressure Drop vs Film Percentage for Variable Diameters")
plt.legend()
plt.show()

########## USE DESMOS TO FIND PRESSURE DROP WITH 30% FILM AT CERTAIN DIAMETER ###########
dp = 17 * PSI2PA # Pressure drop when using 0.3 as smallest diameter
film_percent = 0
film_percent_history = []
m_dot_history = []
diameter_history = []

while film_percent <= 30:
    m_dot = m_dot_ipa * film_percent / 100
    area = m_dot / (cd * np.sqrt(2 * density * dp))
    diameter = np.sqrt(area / (np.pi)) * 2
    m_dot_history.append(m_dot)
    film_percent_history.append(film_percent)
    diameter_history.append(diameter)
    film_percent += 1

plt.plot(film_percent_history, np.array(diameter_history) * M2IN)
plt.xlabel("Film Percent (%)")
plt.ylabel("Diameter (in)")
plt.title("Diameter of Orifice Plate vs Film %")
plt.show()