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
max_film_percent = 15 # [%]
m_dot_ipa = vehicle.parameters.core_fuel_mass_flow_rate # [kg/s]
chamber_pressure = vehicle.parameters.chamber_pressure # [Pa]
manifold_pressure = chamber_pressure / 0.8 # [Pa] 20% Stiffness to prevent backflow
dp = manifold_pressure - chamber_pressure # [Pa]
cd = 0.6
max_m_dot = m_dot_ipa * max_film_percent / 100
density = DENSITY_IPA
area_total = max_m_dot / (cd * np.sqrt(2 * density * dp))
print("STEP ONE: Calculating Total Area of Orifices")
print(f"MAX FILM PERCENT: {max_film_percent:.1f}%")
print(f"TOTAL FILM INJECTION AREA: {area_total * M22IN2:.3f} in^2\n")

film_percent_history_orifices = []
dp_history_orifices = []
for film_percent in np.linspace(0.5, 30, 100):
    m_dot = m_dot_ipa * film_percent / 100
    dp = ((m_dot / (area_total * cd))**2) / (2 * density)
    dp_history_orifices.append(dp)
    film_percent_history_orifices.append(film_percent)
plt.plot(film_percent_history_orifices, np.array(dp_history_orifices) * PA2PSI)
plt.xlabel("Film %")
plt.ylabel("Pressure Drop [psi]")
plt.title(f"Film % vs Pressure Drop\nFixed total area = {area_total * M22IN2:.3f}in^2\nIgnore values, look at shape")
plt.grid()
plt.show()

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

dp_error_percent = 0.1 # percent error from desired pressure drop acceptable
d = 0.8 * IN2M # Diameter of manifold hole 
d0 = 0.00001 * IN2M # Diameter of orifice plate
# di = 0.18 * IN2M # Inner Diameter of tube before metering orifice
# area_i = np.pi * di**2 / 4 # Inner area of tube before metering orifice # NUMBER 1
upstream_pressure = 350 * PSI2PA # Pressure before metering orifice, to be changed
dp = np.inf # Initializing pressure drop
film_percent = 1 # Initializing
film_percent_history = []
m_dot_film_history = []
dp_history = []
d0_history = []
K0_history = []

print("Please wait this might take a second....")

##### Find d0 for all mass flow rates for which desired manifold pressure is achieved #####
target_dp = upstream_pressure - manifold_pressure
dp_error_allowed = (dp_error_percent / 100) * target_dp
for film_percent in np.linspace(3, max_film_percent, 100):
    m_dot_film = m_dot_ipa * film_percent / 100
    #line_velocity = m_dot_film / (density * area_i) # NUMBER 1
    d0 = 0.00001 * IN2M
    dp = np.inf
    while d0 < d:
        K0 = calc_K_sharp_edged_orifice(d, d0)
        A0 = np.pi * d0**2 / 4 # NUMBER 2
        v_orifice = m_dot_film / (density * A0) # NUMBER 2
        dp = K0 * 0.5 * density * v_orifice**2 # NUMBER 2
        #dp = density * (line_velocity**2 / 2) * K0 * (d**4 / d0**4) # NUMBER 1
        err = abs(dp - target_dp)
        if err < dp_error_allowed:
            break

        d0 += 0.00001 * IN2M
    film_percent_history.append(film_percent)
    m_dot_film_history.append(m_dot_film)
    dp_history.append(dp)
    d0_history.append(d0)


plt.plot(np.array(d0_history) * M2IN, film_percent_history)
plt.xlabel("Diameter of Orifice [in]")
plt.ylabel("Film %")
plt.title(f"Diameter of Orifice vs Film %\n[Manifold Diameter = {d * M2IN:.3f}in, Upstream Pressure = {upstream_pressure * PA2PSI:.3f}psi, dp = {dp * PA2PSI:.3f}psi]")
plt.grid()
plt.show()

plt.plot(film_percent_history, np.array(dp_history) * PA2PSI)
plt.xlabel("Film %")
plt.ylabel("Pressure Drop [psi]")
plt.axhline(y = target_dp * PA2PSI, color = 'r')
plt.title("Film % vs Pressure Drop")
plt.grid()
plt.show()

plt.plot(np.array(d0_history) * M2IN, np.array(dp_history) * PA2PSI)
plt.xlabel("Diameter of Orifice [in]")
plt.ylabel("Pressure Drop [psi]")
plt.axhline(y = target_dp * PA2PSI, color = 'r')
plt.title("Diameter of Orifice vs Pressure Drop [psi]")
plt.grid()
plt.show()


##### Validation using cd equation #####
cd = 1/np.sqrt(K0)
testing_history = []
for dp, d0 in zip(dp_history, d0_history):
    A0 = np.pi * d0**2 / 4
    mdot = cd * A0 * np.sqrt(2 * density * dp)
    film_percent = 100 * mdot / m_dot_ipa
    testing_history.append(film_percent)

plt.plot(film_percent_history, testing_history)
plt.title("Ensure that x = y")
plt.xlabel("Film % Calculated")
plt.ylabel("What film % Should be")
plt.grid()
plt.show()

#### testing cd and K stuff ####
d = 0.8 * IN2M
d0 = 0.01 * IN2M
d0_history_test = []
K0_history_test = []
cd_history_test = []
while d0 <= d:
    K0 = calc_K_sharp_edged_orifice(d, d0)
    cd = 1/np.sqrt(K0)
    d0_history_test.append(d0)
    K0_history_test.append(K0)
    cd_history_test.append(cd)
    d0 += 0.01 * IN2M

plt.plot(np.array(d0_history_test) * M2IN, K0_history_test)
plt.xlabel("Orifice Diameter (in)")
plt.ylabel("K")
plt.title(f"Manifold Diameter: {d * M2IN}in")
plt.show()

plt.plot(np.array(d0_history_test) * M2IN, cd_history_test)
plt.xlabel("Orifice Diameter (in)")
plt.ylabel("Cd")
plt.title(f"Manifold Diameter: {d * M2IN}in")
plt.show()

plt.plot(K0_history_test, cd_history_test)
plt.xlabel("K")
plt.ylabel("Cd")
plt.title(f"Manifold Diameter: {d * M2IN}in")
plt.show()