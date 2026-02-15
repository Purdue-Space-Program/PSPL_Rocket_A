# Stolen from Copperhead

# Press Cv Calcs
# This code conservatively estimates the maximum required Cv for our pressurization system (applicable to both regulators and solenoid valves)
# Owned by Hugo Filmer

from CoolProp.CoolProp import PropsSI
from math import sqrt

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
from vehicle_parameters import parameters


def Cv_Choked(press_gas, m_dot, P1, T1):

    # Unit conversions
    c.PSI2PA = 6894.76
    c.PA2PSI = 1 / c.PSI2PA
    c.M32FT3 = 35.3147
    SEC2MIN = 1 / 60
    K2R = 9 / 5

    # Constants
    N2 = 22.67 # units constant for flow rate equation

    # Standard conditions
    P_STD = 14.7 * c.PSI2PA # [pa] standard pressure
    T_STD = 288.7 # [K] standard temperature

    rho_std_air = PropsSI('D', 'P', P_STD, 'T', T_STD, 'air') # [kg/^3] standard density of air
    rho_std_press_gas = PropsSI('D', 'P', P_STD, 'T', T_STD, press_gas) # [kg/m^3] standard density of press gas
    
    Gg_press_gas = rho_std_press_gas / rho_std_air # [] specific gravity of press gas

    q_std_press_gas = m_dot / rho_std_press_gas # [m^3/s] standard volumetric flow rate of press gas
    scfm_press_gas = q_std_press_gas * c.M32FT3 / SEC2MIN # [ft^3/min] standard volumetric flow rate of press gas, in SCFM

    P1_psi = P1 * c.PA2PSI # [psia] inlet pressure in psia
    T1_r = T1 * K2R # [deg R] inlet temperature in deg R

    cv = scfm_press_gas * sqrt(Gg_press_gas * T1_r) / (0.471 * N2 * P1_psi)

    return cv


# Unit conversions
L2M3 = 0.001
PSI2PA = 6894.76
LBM2KG = 0.453592

# Parameters
# Press system
V_COPV = parameters.COPV_volume # [m^3] COPV volume
T1_COPV = 300 # [K] starting COPV temperature (assumed)
PRESS_GAS = 'nitrogen'
# Propellant system
# Oxidizer
P_FILL = 40 * PSI2PA # [Pa] fill pressure for LOx (assumed)
OXIDIZER_NAME = 'oxygen'
# Fuel
M_DOT_FU = parameters.core_fuel_mass_flow_rate + parameters.film_fuel_mass_flow_rate # [kg/s] fuel mass flow rate
FUEL_NAME = 'ethanol'

# Volumetric flow rates
T_ox = PropsSI('T', 'P', P_FILL, 'Q', 0, OXIDIZER_NAME) # [K] oxidizer temeprature
rho_ox = PropsSI('D', 'P', parameters.oxidizer_tank_pressure, 'T', T_ox, OXIDIZER_NAME) # [kg/m^3] oxidizer density
Q_dot_ox = parameters.oxidizer_mass_flow_rate / rho_ox # [m^3/s] oxidizer volumetric flow rate

rho_fu = PropsSI('D', 'P', parameters.fuel_tank_pressure, 'T', c.T_AMBIENT, FUEL_NAME) # [kg/m^3] fuel density
Q_dot_fu = parameters.oxidizer_mass_flow_rate / rho_fu # [m^3/s] fuel volumetric flow rate

Q_dot_ox_press = Q_dot_ox # [m^3/s] volumetric flow rate of pressurant required to pressurize the oxidizer tank
Q_dot_fu_press = Q_dot_fu # [m^3/s] volumetric flow rate of pressurant required to pressurize the fuel tank

# Mass flow rates
# Assume worst-case total collapse for press gas in ox tank
rho_ox_press = PropsSI('D', 'P', parameters.oxidizer_tank_pressure, 'T', T_ox, PRESS_GAS) # [kg/m^3] assumed density of press gas in ox tank
m_dot_ox_press = rho_ox_press * Q_dot_ox_press # [kg/s] required mass flow rate of press gas into ox tank
# Assume end-of-burn conditions with isentropic expansion in COPV and adiabatic upper plumbing for coldest press gas going into fuel tank

P2_copv = 1000 * c.PSI2PA # (2 * max(P_OX, P_FU)) # [Pa] ending COPV pressure (assumed, should vary this to check numbers)
s1_copv = PropsSI('S', 'P', parameters.COPV_starting_pressure, 'T', T1_COPV, PRESS_GAS) # [J/kgK] starting COPV entropy
s2_copv = s1_copv # [J/kgK] ending COPV entropy (assumed isentropic expansion)
h2_copv = PropsSI('H', 'P', P2_copv, 'S', s2_copv, PRESS_GAS) # [J/kg] ending COPV enthalpy
h2_fu_press = h2_copv # [J/kg] ending fuel tank inlet enthalpy (stagnation enthalpy is constant for an adiabatic flow with no work)
rho_fu_press = PropsSI('D', 'P', parameters.fuel_tank_pressure, 'H', h2_fu_press, PRESS_GAS) # [kg/m^3] assumed density of press gas in fuel tank
m_dot_fu_press = rho_fu_press * Q_dot_fu_press # [kg/s] required mass flow rate of press gas into fuel tank

# Required Cv calcuations
T2_copv = PropsSI('T', 'P', P2_copv, 'S', s2_copv, PRESS_GAS) # [K] ending COPV temperature
# Ox tank
cv_required_ox = Cv_Choked(PRESS_GAS, m_dot_ox_press, P2_copv, T2_copv)
# Fuel tank
cv_required_fu = Cv_Choked(PRESS_GAS, m_dot_fu_press, P2_copv, T2_copv)

# Print output
print(f"m_dot_ox_press: {m_dot_ox_press:.2f}")
print(f"m_dot_fu_press: {m_dot_fu_press:.2f}")

print(f"The maximum required Cv for pressurizing the oxidizer tank is: {cv_required_ox:.5f}")
print(f"The maximum required Cv for pressurizing the fuel tank is:     {cv_required_fu:.5f}")