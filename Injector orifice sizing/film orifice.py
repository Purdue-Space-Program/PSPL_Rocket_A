from CoolProp.CoolProp import PropsSI
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *
import vehicle_parameters as vehicle

def find_film_mdot(m_dot_ipa, film_percent):
    m_dot = film_percent * m_dot_ipa / (1 - (film_percent / 100))
    return m_dot
##### Find total area of the orifices needed to sustain a maximum of 30% film #####
max_film_percent = 30 # [%]
m_dot_ipa = vehicle.parameters.core_fuel_mass_flow_rate # [kg/s]
p1 = 310 * PSI2PA # [Pa]
p2 = vehicle.parameters.chamber_pressure # [Pa]
dp = p1 - p2 # [Pa]
cd = 0.6
max_m_dot = find_film_mdot(m_dot_ipa, max_film_percent)
density = DENSITY_IPA
area_total = max_m_dot / (cd * np.sqrt(2 * density * dp))
print("STEP ONE: Calculating Total Area of Orifices")
print(f"MAX FILM PERCENT: {max_film_percent * 100:.1f}%")
print(f"TOTAL AREA: {area_total * M22IN2:.3f} in^2\n")

##### Graph dp through adjustable tap-off orifice and film mass flow rate #####
film_percent = 1
film_percent_history = [0, ]
pressure_drop_history = [0, ]
m_dot_history = [0, ]
while film_percent <= max_film_percent:
    m_dot = find_film_mdot(m_dot_ipa, max_film_percent)
    m_dot_history.append(m_dot)
    film_percent_history.append(film_percent)
    dp = 0 
    pressure_drop_history.append(dp)
    film_percent += 1

print("STEP 2: Graphing Pressure Drop Through Orifice with Film Mass Flow Rate")