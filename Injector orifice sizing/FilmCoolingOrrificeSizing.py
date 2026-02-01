import numpy as np
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import vehicle_parameters as vp
import matplotlib.pyplot as plt
import constants as c

m_dotTotal = vp.parameters.total_mass_flow_rate #total propellant mass flow
f_film = 0.20 #film cooling fraction
m_dot_film = m_dotTotal * f_film #mass flow rate for film cooling

#Constants of film gas
gamma_film = 10 #film cooling specific heat ratio
M = 1 #molar mass of film cooling gas (kg/mol)
R = 8.314 / M #specific gas constant for film cooling gas (J/kg-K)
T0 = 300 #stagnation temperature of film cooling gas (K)
C_d = 0.75 #discharge coefficient of orifice

P_main = vp.parameters.chamber_pressure #combustion chamber pressure (Pa)
P_film = 2 #film manifold pressure (Pa)

#check for choked flow

pressure_ratio_critical = P_main / P_film
critical_ratio = (2 / (gamma_film + 1)) ** (gamma_film / (gamma_film - 1))

if pressure_ratio_critical > critical_ratio: #not choked
    term_1 = (pressure_ratio_critical ** (2/gamma_film)) - (pressure_ratio_critical ** ((gamma_film + 1)/gamma_film))
    term_2 = (2 * gamma_film) / (R * T0 * (gamma_film - 1))
    Area_orifice = m_dot_film / (C_d * P_main * np.sqrt(term_1 * term_2))

if pressure_ratio_critical <= critical_ratio: #choked flow
    term_1 = (2 / (gamma_film + 1)) ** ((gamma_film + 1) / (gamma_film - 1))
    term_2 = gamma_film / (R * T0)
    Area_orifice = m_dot_film / (C_d * P_main * np.sqrt(term_1 * term_2))

diameter_orifice = 2 * np.sqrt(Area_orifice / np.pi)