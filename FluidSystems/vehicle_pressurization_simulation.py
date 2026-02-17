### This code was adapted from the sim made by Hugo Filmer to work for Rocket A. David Gustafsson is the one who adapted this

# This code simulates the pressurization process in the propellant tanks over time, assuming a constant pressure in each tank.
# Its primary use is to determine the final COPV pressure and therefore if the COPV is large enough.
# Hugo Filmer

from CoolProp.CoolProp import PropsSI
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from heat_transfer_functions import *

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
from vehicle_parameters import parameters as p

# Simulation settings
T_AMBIENT = 293 # [K] ambient temperature
LOITER_TIME = 5 * 1 # [s] time between prepressurization and the start of flow
LAG_TIME = 0 # [s] time the simulation should continue to run for after the run valves are closed
DT = 0.050 # [s] simulation step size
Q_FACTOR = 1 # [] factor to multiply heat transfer by (for conservatism)
TEXT_OUTPUT = True # True to print summary output, including conservation and EoS checks
PLOT_OUTPUT = True # True to make pretty plots of the results :)

# inputs

# Vehicle parameters
# COPV
vehicle_name = "Rocket_A"

if vehicle_name == "Copperhead":
    ADIABATIC = False # True to ignore heat transfer in the tanks
    PREPRESS = "isothermal" # prepressurization type (choose from 'adiabatic' or 'isothermal') isothermal means infinite loiter time
    PRESS_LINE_CHILL = False # True to account for heat transfer in the helium line that runs through the oxidizer tank
    
    GRAVITY = 1 * 9.81 # [m/s/s] local gravitational acceleration (may be > 9.81 in flight)
    # COPV
    PRESS_GAS = 'helium'
    P_COPV = 4500 * c.PSI2PA # [Pa] starting COPV pressure
    T_COPV = 300 # [K] starting COPV temperature (assumed)
    V_COPV = 12 * c.L2M3 # [m^3] COPV volume
    # Tanks
    D_TANK = 8.84 * c.IN2M # [m] tank outer diameter
    T_TANK = 0.17 * c.IN2M # [m] tank wall thickness
    RHO_TANK = 2700 # [kg/m^3] tank material density (aluminum)
    M_BULKHEAD = 1.5 * c.LBM2KG # [kg] bulkhead mass (single bulkhead)
    CP_TANK = 500 # [J/kgK] tank material specific heat
    D_PRESS_LINE = 3/8 * c.IN2M # [m] fuel tank pressurization line outer diameter
    T_PRESS_LINE = 0.049 * c.IN2M # [m] fuel tank pressurization line wall thickness
    P_TANK = 287 * c.PSI2PA
    # Oxidizer
    OXIDIZER = 'oxygen'
    P_FILL = 40 * c.PSI2PA # [Pa] fill pressure for LOx (assumed)
    M_DOT_OX = 5.503 * c.LBM2KG # [kg/s] LOx mass flow rate
    P_OX = P_TANK # [Pa] oxidizer tank nominal pressure
    V_OX = 2246 * c.IN32M3 # [m^3] oxidizer tank total volume
    ULLAGE_OX = 1 / 100 # [] oxidizer tank volume ullage fraction
    RESIDUAL_OX = 1 / 100 # [] oxidizer tank volume residual fraction
    # note: drain time is based on oxidizer residuals since that's what we'll be sensing
    # Fuel
    FUEL = 'ethanol'
    M_DOT_FU = (3.336 + 0.500) * c.LBM2KG # [kg/s] fuel mass flow rate
    P_FU = P_TANK # [Pa] fuel tank nominal pressure
    V_FU = 2062 * c.IN32M3 # [m^3] fuel tank total volume
    ULLAGE_FU =  1 / 100 # [] fuel tank volume ullage fraction
    
    V_ullage_ox = 98.9 * c.IN32M3 # [m^3] oxidizer tank ullage volume
    V_ullage_fu = 10 * c.IN32M3 # [m^3] fuel tank ullage volume
    
    # V_ullage_ox = V_OX * ULLAGE_OX # [m^3] oxidizer tank ullage volume
    # V_ullage_fu = V_FU * ULLAGE_FU # [m^3] fuel tank ullage volume

if vehicle_name == "Rocket_A":
    ADIABATIC = False # True to ignore heat transfer in the tanks
    PREPRESS = "isothermal" # prepressurization type (choose from 'adiabatic' or 'isothermal') isothermal means infinite loiter time
    PRESS_LINE_CHILL = False # True to account for heat transfer in the helium line that runs through the oxidizer tank
    
    PRESS_GAS = "nitrogen"
    P_COPV = p.COPV_starting_pressure # [Pa] starting COPV pressure
    T_COPV = 300 # [K] starting COPV temperature (assumed)
    GRAVITY = p.one_DoF_off_the_rail_acceleration * 9.81 # [m/s/s] local gravitational acceleration (may be > 9.81 in flight)
    V_COPV = p.COPV_volume # [m^3] COPV volume

    # Tanks
    D_TANK = p.tank_outer_diameter # [m] tank outer diameter
    T_TANK = p.tank_wall_thickness # [m] tank wall thickness
    RHO_TANK = c.DENSITY_AL # [kg/m^3] tank material density (aluminum)
    M_BULKHEAD = 1.5 * c.LBM2KG # [kg] bulkhead mass (single bulkhead)
    CP_TANK = 500 # [J/kgK] tank material specific heat
    D_PRESS_LINE = (3/8) * c.IN2M # [m] fuel tank pressurization line outer diameter
    T_PRESS_LINE = 0.049 * c.IN2M # [m] fuel tank pressurization line wall thickness
    P_TANK = p.oxidizer_tank_pressure
    
    # Oxidizer
    OXIDIZER = "oxygen"
    P_FILL = 40 * c.PSI2PA # [Pa] fill pressure for LOx (assumed)
    M_DOT_OX = p.oxidizer_mass_flow_rate # [kg/s] LOx mass flow rate
    P_OX = P_TANK # [Pa] oxidizer tank nominal pressure
    V_OX = p.oxidizer_tank_usable_volume # [m^3] oxidizer tank total volume ############################################################ FIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIX
    ############################################################ FIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIX
    ############################################################ FIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIX
    ############################################################ FIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIX
    ############################################################ FIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIX
    ############################################################ FIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIX
    ############################################################ FIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIX
    ULLAGE_OX = 10 / 100 # [] oxidizer tank volume ullage fraction
    RESIDUAL_OX = 10 / 100 # [] oxidizer tank volume residual fraction
    # note: drain time is based on oxidizer residuals since that's what we'll be sensing
    
    # Fuel
    FUEL = "ethanol"
    M_DOT_FU = p.core_fuel_mass_flow_rate # [kg/s] fuel mass flow rate
    P_FU = P_TANK # [Pa] fuel tank nominal pressure
    V_FU = p.fuel_tank_usable_volume # [m^3] fuel tank total volume ############################################################ FIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIX
    ############################################################ FIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIX
    ############################################################ FIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIX
    ############################################################ FIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIX
    ############################################################ FIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIX
    ############################################################ FIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIX
    ############################################################ FIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIXFIX
    ULLAGE_FU =  3 / 100 # [] fuel tank volume ullage fraction
    
    V_ullage_ox = V_OX * ULLAGE_OX # [m^3] oxidizer tank ullage volume
    V_ullage_fu = V_FU * ULLAGE_FU # [m^3] fuel tank ullage volume
else:
    raise ValueError("what")

# Volumetric flow rates
# Oxidizer
T_fill_ox = PropsSI('T', 'P', P_FILL, 'Q', 0, OXIDIZER) # [K] oxidizer temperature
rho_ox_nom = PropsSI('D', 'P', P_OX, 'T', T_fill_ox, OXIDIZER) # [kg/m^3] oxidizer density
V_dot_ox_nom = M_DOT_OX / rho_ox_nom # [m^3/s] oxidizer volumetric flow rate
# Fuel
T_fill_fu = T_AMBIENT # [K] fuel temperature
rho_fu_nom = PropsSI('D', 'P', P_FU, 'T', T_AMBIENT, FUEL) # [kg/m^3] fuel density
V_dot_fu_nom = M_DOT_FU / rho_fu_nom # [m^3/s] fuel volumetric flow rate

# Gas properties
R_press_gas = 8.314462 / PropsSI('M', PRESS_GAS) # [J/kgK] pressurizing gas constant

# Tank properties
D_tank_inner = D_TANK - (2 * T_TANK)
A_bulkhead = np.pi * (D_tank_inner**2) / 4 # [m^2] bulkhead area (flat circle)
l_tank_ox = V_OX / A_bulkhead # [m] oxidizer tank length
m_tank_ox = 2 * M_BULKHEAD + np.pi / 4 * (D_TANK**2 - D_tank_inner**2) * l_tank_ox * RHO_TANK # [kg] oxidizer tank mass
l_tank_fu = V_FU / A_bulkhead # [m] fuel tank length
m_tank_fu = 2 * M_BULKHEAD + np.pi / 4 * (D_TANK**2 - D_tank_inner**2) * l_tank_fu * RHO_TANK # [kg] fuel tank mass
D_inner_press_line = D_PRESS_LINE - 2 * T_PRESS_LINE # [m] fuel tank press line inner diameter
A_inner_press_line = np.pi * (D_inner_press_line**2) / 4 # [m^2] fuel tank press line inner area
L_press_line = l_tank_ox # [m] fuel tank press line length



# Starting COPV pressure and temperature, as well as oxidizer and fuel tank temperature and number of iterations needed, after  pre-pressurization
if PREPRESS == "adiabatic": # Calculate initial conditions assuming adiabatic prepress
    [P0_copv, T0_copv, T0_ox, T0_fu, iter_count] = adiabatic_press(P_COPV, T_COPV, V_COPV, P_OX, V_ullage_ox, P_FU, V_ullage_fu, PRESS_GAS)
elif PREPRESS == "isothermal": # Calculate initial conditions assuming tank ullage temperatures equal tank temperatures (infinite loiter time)
    [P0_copv, T0_copv, T0_ox, T0_fu, iter_count] = isothermal_press(P_COPV, T_COPV, V_COPV, P_OX, V_ullage_ox, T_fill_ox, P_FU, V_ullage_fu, T_fill_fu, PRESS_GAS)
else:
    raise Exception("Invalid prepressurization model defined. Please choose from 'adiabatic' or 'isothermal'.")

s_copv = PropsSI('S', 'P', P0_copv, 'T', T0_copv, PRESS_GAS) # [J/kgK] COPV specific entropy (constant)
e0_ox = PropsSI('U', 'P', P_OX, 'T', T0_ox, PRESS_GAS) # [J/kgK] oxidizer tank specific energy
e0_fu = PropsSI('U', 'P', P_FU, 'T', T0_fu, PRESS_GAS) # [J/kgK] fuel tank specific energy

if TEXT_OUTPUT == True:
    print('INITIAL CONDITIONS -------------------------------------------------------------------------------------------')
    print(f'The starting COPV conditions are {P0_copv/c.PSI2PA:.3f} [psi] and {T0_copv:.3f} [K].')
    print(f'The starting tank temperatures are {T0_ox:.3f} [K] in ox and {T0_fu:.3f} [K] in fuel.')

# Simulation time
V_residual_ox = V_OX * RESIDUAL_OX # [m^3] oxidizer tank residual volume

drain_time = (V_OX - V_ullage_ox - V_residual_ox) / V_dot_ox_nom # [s] time to drain the oxidizer tank to the residual volume
# drain_time = parameters.burn_time # [s] time to drain the oxidizer tank to the residual volume

total_time = LOITER_TIME + drain_time + LAG_TIME # [s] total time the simulation will run for
total_steps = int(total_time // DT) # [] total number of simulation steps
times = np.linspace(0, total_time, total_steps + 1) # [s] simulation time array

# Initial conditions
# COPV
P_copv = np.empty(total_steps + 1) # [Pa] COPV pressure array
P_copv[0] = P0_copv
T_copv = np.empty(total_steps + 1) # [K] COPV temperature array
T_copv[0] = T0_copv
rho_copv = np.empty(total_steps + 1) # [kg/m^3] COPV density array
rho_copv[0] = PropsSI('D', 'P', P_copv[0], 'T', T_copv[0], PRESS_GAS)
m_copv = np.empty(total_steps + 1) # [kg] COPV mass array
m_copv[0] = V_COPV * rho_copv[0]
# Tank
V_ullage_tank_ox = V_ullage_ox # [m^3] oxidizer tank ullage volume
V_ullage_tank_fu = V_ullage_fu # [m^3] fuel tank ullage volume
T_tank_wall_ox = np.empty(total_steps + 1) # [K] oxidizer tank wall temperature array
T_tank_wall_ox[0] = T_fill_ox
T_tank_wall_fu = np.empty(total_steps + 1) # [K] fuel tank wall temperature array
T_tank_wall_fu[0] = T_AMBIENT
# Heat transfer coefficients
# In oxidizer tank
h_ox_gas_wall = np.empty(total_steps + 1) # [W/m^2K] convective heat transfer coefficient between oxidizer ullage gas and tank wall array
h_ox_liquid_wall = np.empty(total_steps + 1) # [W/m^2K] convective heat transfer coefficient between oxidizer propellant and tank wall array
h_ox_gas_plate = np.empty(total_steps + 1) # [W/m^2K] convective heat transfer coefficient between oxidizer ullage gas and upper bulkhead array
h_ox_liquid_plate = np.empty(total_steps + 1) # [W/m^2K] convective heat transfer coefficient between oxidizer propellant and lower bulkhead array
h_ox_interface = np.empty(total_steps + 1) # [W/m^2K] convective heat transfer coefficient between oxidizer ullage gas and oxidizer propellant array
# In fuel tank
h_fu_gas_wall = np.empty(total_steps + 1) # [W/m^2K] convective heat transfer coefficient between fuel ullage gas and tank wall array
h_fu_liquid_wall = np.empty(total_steps + 1) # [W/m^2K] convective heat transfer coefficient between fuel propellant and tank wall array
h_fu_gas_plate = np.empty(total_steps + 1) # [W/m^2K] convective heat transfer coefficient between fuel ullage gas and upper bulkhead array
h_fu_liquid_plate = np.empty(total_steps + 1) # [W/m^2K] convective heat transfer coefficient between fuel propellant and lower bulkhead array
h_fu_interface = np.empty(total_steps + 1) # [W/m^2K] convective heat transfer coefficient between fuel ullage gas and fuel propellant array
# Averages
h_avg_ullage_ox = np.empty(total_steps + 1) # [W/m^2K] average convective heat transfer coefficient into oxidizer ullage gas for plotting
h_avg_ullage_fu = np.empty(total_steps + 1) # [W/m^2K] average convective heat transfer coefficient into fuel ullage gas for plotting
# Oxidizer
e_ullage_ox = np.empty(total_steps + 1) # [J/kg] oxidizer ullage energy array
e_ullage_ox[0] = e0_ox
T_ullage_ox = np.empty(total_steps + 1) # [K] oxidizer ullage temperature array
T_ullage_ox[0] = T0_ox
m_ullage_ox = np.empty(total_steps + 1) # [kg] oxidizer ullage mass array
m_ullage_ox[0] = V_ullage_tank_ox * PropsSI('D', 'P', P_OX, 'T', T0_ox, PRESS_GAS)
mdot_ullage_ox = np.empty(total_steps + 1) # [kg/s] oxidizer ullage mass flow rate array
rho_ullage_ox = np.empty(total_steps + 1) # [kg/m^3] oxidizer ullage density array
rho_ullage_ox[0] = m_ullage_ox[0] / V_ullage_tank_ox
partial_derivative_of_ullage_density_with_respect_to_internal_energy_in_ox = np.empty(total_steps + 1) # [1/Jm^3] partial derivative of oxidizer ullage density WRT internal energy array
partial_derivative_of_ullage_density_with_respect_to_internal_energy_in_ox[0] = PropsSI('d(D)/d(U)|P', 'D', rho_ullage_ox[0], 'U', e_ullage_ox[0], PRESS_GAS) # [1/Jm^3] partial derivative of oxidizer ullage density WRT internal energy


# Fuel
e_ullage_fu = np.empty(total_steps + 1) # [J/kg] fuel ullage energy array
e_ullage_fu[0] = e0_fu
T_ullage_fu = np.empty(total_steps + 1) # [K] fuel ullage temperature array
T_ullage_fu[0] = T0_fu
m_ullage_fu = np.empty(total_steps + 1) # [kg] fuel ullage mass array
m_ullage_fu[0] = V_ullage_tank_fu * PropsSI('D', 'P', P_FU, 'T', T0_fu, PRESS_GAS)
mdot_ullage_fu = np.empty(total_steps + 1) # [kg/s] fuel ullage mass flow rate array
rho_ullage_fu = np.empty(total_steps + 1) # [kg/m^3] fuel ullage density array
rho_ullage_fu[0] = m_ullage_fu[0] / V_ullage_tank_fu

# Press line
h_press_line = np.empty(total_steps + 1) # [W/m^2K] fuel pressurization line convective heat transfer coefficient array
Q_dot_press_line = np.empty(total_steps + 1) # [J/s] fuel pressurization line heat transfer rate array

# Simulation loop
step = 0 # current simulation step
time = 0 # [s] current simulation time

for i in tqdm(range(total_steps)):

    # Calculate heat transfers
    if ADIABATIC == True: # No heat transfer :)
        Q_dot_ox = 0
        Q_dot_fu_ullage = 0
        Q_dot_ox_tank = 0
        Q_dot_fu_tank = 0
    else:
        # Oxidizer coefficients
        L_char_ox = V_ullage_tank_ox / A_bulkhead # [m] Characteristic length for tank wall
        h_ox_gas_wall[step] = wall_convect(GRAVITY, P_OX, T_tank_wall_ox[step], T_ullage_ox[step], L_char_ox, PRESS_GAS) # [W/m^2K] heat transfer coefficient between oxidizer ullage and tank wall
        h_ox_liquid_wall[step] = wall_convect(GRAVITY, P_OX, T_tank_wall_ox[step], T_fill_ox, L_char_ox, OXIDIZER) # [W/m^2K] heat transfer coefficient between oxidizer and tank wall
        h_ox_gas_plate[step] = plate_convect(GRAVITY, P_OX, T_tank_wall_ox[step], T_ullage_ox[step], D_tank_inner/4, PRESS_GAS, 'enhanced') # [W/m^2K] heat transfer coefficient between oxidizer ullage and tank bulkhead
        h_ox_liquid_plate[step] = plate_convect(GRAVITY, P_OX, T_tank_wall_ox[step], T_fill_ox, D_tank_inner/4, OXIDIZER, 'enhanced') # [W/m^2K] heat transfer coefficient between oxidizer and tank bulkhead
        h_ox_interface[step] = interface_convect(GRAVITY, P_OX, T_ullage_ox[step], T_fill_ox, D_tank_inner/4, PRESS_GAS, OXIDIZER) # [W/m^2K] heat transfer coefficient between oxidizer ullage and oxidizer
        # Oxidizer heat transfers
        # Areas
        A_ox_gas_wall = np.pi*D_tank_inner*L_char_ox # [m^2] oxidizer tank wall area exposed to ullage gas  
        A_ox_liquid_wall = np.pi*D_tank_inner*(V_OX - V_ullage_tank_ox) / A_bulkhead # [m^2] oxidizer tank wall area exposed to oxidizer 
        # Q_dots
        Q_dot_ox_gas_wall = h_ox_gas_wall[step] * (A_ox_gas_wall) * (T_tank_wall_ox[step] - T_ullage_ox[step]) # [W] heat transfer rate between oxidizer ullage and tank wall
        Q_dot_ox_liquid_wall = h_ox_liquid_wall[step] * (A_ox_liquid_wall) * (T_fill_ox - T_tank_wall_ox[step]) # [W/m^2K] heat transfer rate between oxidizer and tank wall
        Q_dot_ox_gas_plate = h_ox_gas_plate[step]  * (A_bulkhead) * (T_tank_wall_ox[step] - T_ullage_ox[step]) # [W/m^2K] heat transfer rate between oxidizer ullage and tank bulkhead
        Q_dot_ox_liquid_plate = h_ox_liquid_plate[step] * (A_bulkhead) * (T_fill_ox - T_tank_wall_ox[step]) # [W/m^2K] heat transfer rate between oxidizer and tank bulkhead
        Q_dot_ox_interface = h_ox_interface[step] * (A_bulkhead) * (T_fill_ox - T_ullage_ox[step]) # [W/m^2K] heat transfer rate between oxidizer ullage and oxidizer

        Q_dot_ox = Q_FACTOR * (Q_dot_ox_gas_wall + Q_dot_ox_gas_plate + Q_dot_ox_interface) # [W] total heat transfer from oxidizer ullage
        Q_dot_ox_tank = Q_FACTOR * (Q_dot_ox_liquid_wall + Q_dot_ox_liquid_plate) # [W] total heat transfer to oxidizer tank

        # Fuel coefficients
        L_char_fu = V_ullage_tank_fu / A_bulkhead # [m] Characteristic length for tank wall
        h_fu_gas_wall[step] = wall_convect(GRAVITY, P_FU, T_tank_wall_fu[step], T_ullage_fu[step], L_char_fu, PRESS_GAS)  # [W/m^2K] heat transfer coefficient between fuel ullage and tank wall
        h_fu_liquid_wall[step] = wall_convect(GRAVITY, P_FU, T_tank_wall_fu[step], T_fill_fu, L_char_fu, FUEL) # [W/m^2K] heat transfer coefficient between fuel and tank wall
        if T_ullage_fu[step] > T_tank_wall_fu[step]:
            enhancement = 'enhanced'
        else:
            enhancement = 'reduced'
        h_fu_gas_plate[step] = plate_convect(GRAVITY, P_FU, T_tank_wall_fu[step], T_ullage_fu[step], D_tank_inner/4, PRESS_GAS, enhancement) # [W/m^2K] heat transfer coefficient between fuel ullage and tank bulkhead
        h_fu_liquid_plate[step] = plate_convect(GRAVITY, P_FU, T_tank_wall_fu[step], T_fill_fu, D_tank_inner/4, FUEL, enhancement) # [W/m^2K] heat transfer coefficient between fuel and tank bulkhead
        h_fu_interface[step] = interface_convect(GRAVITY, P_FU, T_ullage_fu[step], T_fill_fu, D_tank_inner/4, PRESS_GAS, FUEL) # [W/m^2K] heat transfer coefficient between fuel ullage and oxidizer
        # Fuel heat transfers
        # Areas
        A_fu_gas_wall = np.pi*D_tank_inner*L_char_fu # [m^2] fuel tank wall area exposed to ullage gas
        A_fu_liquid_wall = np.pi*D_tank_inner*(V_FU - V_ullage_tank_fu) / A_bulkhead # [m^2] fuel tank wall area exposed to oxidizer
        # Q_dots
        Q_dot_fu_gas_wall = h_fu_gas_wall[step] * (A_fu_gas_wall) * (T_tank_wall_fu[step] - T_ullage_fu[step]) # [W] heat transfer rate between fuel ullage and tank wall
        Q_dot_fu_liquid_wall = h_fu_liquid_wall[step] * (A_fu_liquid_wall) * (T_fill_fu - T_tank_wall_fu[step]) # [W/m^2K] heat transfer rate between fuel and tank wall
        Q_dot_fu_gas_plate = h_fu_gas_plate[step]  * (A_bulkhead) * (T_tank_wall_fu[step] - T_ullage_fu[step]) # [W/m^2K] heat transfer rate between fuel ullage and tank bulkhead
        Q_dot_fu_liquid_plate = h_fu_liquid_plate[step] * (A_bulkhead) * (T_fill_fu - T_tank_wall_fu[step])  # [W/m^2K] heat transfer rate between fuel and tank bulkhead
        Q_dot_fu_interface = h_fu_interface[step] * (A_bulkhead) * (T_fill_fu - T_ullage_fu[step]) # [W/m^2K] heat transfer rate between fuel ullage and oxidizer

        Q_dot_fu_ullage = Q_FACTOR * (Q_dot_fu_gas_wall + Q_dot_fu_gas_plate + Q_dot_fu_interface) # [W] total heat transfer from fuel ullage
        Q_dot_fu_tank = Q_FACTOR * (Q_dot_fu_liquid_wall + Q_dot_fu_liquid_plate) # [W] total heat transfer to fuel tank

        # Average coefficients for plotting
        h_avg_ullage_ox[step] = Q_dot_ox / ((2*A_bulkhead + A_ox_gas_wall) * (T_tank_wall_ox[step] - T_ullage_ox[step])) # [W/m^2K] averaged oxidizer ullage heat transfer coefficient
        h_avg_ullage_fu[step] = Q_dot_fu_ullage / ((2*A_bulkhead + A_fu_gas_wall) * (T_tank_wall_fu[step] - T_ullage_fu[step])) # [W/m^2K] averaged fuel ullage heat transfer coefficient

    # Outflow rate
    if LOITER_TIME <= time < LOITER_TIME + drain_time:
        V_dot_ox = V_dot_ox_nom # [m^3/s] volumetric flow rate of propellant out of oxidizer tank
        V_dot_fu = V_dot_fu_nom # [m^3/s] volumetric flow rate of propellant out of fuel tank
    else:
        V_dot_ox = 0
        V_dot_fu = 0

    # Calculate tank properties
    rho_ullage_ox[step] = m_ullage_ox[step] / V_ullage_tank_ox
    rho_ullage_fu[step] = m_ullage_fu[step] / V_ullage_tank_fu # [kg/m^3] fuel ullage density
    
    h_in = PropsSI('H', 'D', rho_copv[step], 'S', s_copv, PRESS_GAS) # [J/kg] adiabatic flow between COPV and tank so isenthalpic
    try:
        # print(f"\nrho_ullage_ox[step]: {rho_ullage_ox[step]}")
        # print(f"e_ullage_ox[step]: {e_ullage_ox[step]}")
        partial_derivative_of_ullage_density_with_respect_to_internal_energy_in_ox[step] = PropsSI('d(D)/d(U)|P', 'D', rho_ullage_ox[step], 'U', e_ullage_ox[step], PRESS_GAS) # [1/Jm^3] partial derivative of oxidizer ullage density WRT internal energy
        # print(f"partial_derivative_of_ullage_density_with_respect_to_internal_energy_in_ox: {partial_derivative_of_ullage_density_with_respect_to_internal_energy_in_ox[step]}")
    except:
        plt.plot(time, partial_derivative_of_ullage_density_with_respect_to_internal_energy_in_ox)
        plt.show
        
    
    
    # Oxidizer e_dot
    e_dot_ox = (
        (Q_dot_ox - P_OX*V_dot_ox - rho_ullage_ox[step]*V_dot_ox*e_ullage_ox[step] + rho_ullage_ox[step]*V_dot_ox*h_in)
        / (m_ullage_ox[step] + partial_derivative_of_ullage_density_with_respect_to_internal_energy_in_ox[step]*V_ullage_tank_ox*e_ullage_ox[step] - partial_derivative_of_ullage_density_with_respect_to_internal_energy_in_ox[step]*V_ullage_tank_ox*h_in)
    ) # [J/kg/s] rate of change of internal energy in oxidizer ullage
    # Oxidizer Mdot
    mdot_ullage_ox[step] = partial_derivative_of_ullage_density_with_respect_to_internal_energy_in_ox[step]*e_dot_ox*V_ullage_tank_ox + rho_ullage_ox[step]*V_dot_ox # [kg/s] rate of change of mass in oxidizer ullage

    partial_derivative_of_ullage_density_with_respect_to_internal_energy_in_fuel = PropsSI('d(D)/d(U)|P', 'D', rho_ullage_fu[step], 'U', e_ullage_fu[step], PRESS_GAS) # [1/Jm^3] partial derivative of fuel ullage density WRT internal energy

    # Fuel e_dot
    if PRESS_LINE_CHILL == True:
        if step > 0:
            [h_press_line[step], Q_dot_press_line[step]] = line_heat_transfer(P_FU, h_in, T_fill_ox, D_inner_press_line, L_press_line, mdot_ullage_fu[step - 1], PRESS_GAS)
            h_in += (Q_dot_press_line[step] / mdot_ullage_fu[step - 1])

    e_dot_fu = (
        (Q_dot_fu_ullage - P_FU*V_dot_fu - rho_ullage_fu[step]*V_dot_fu*e_ullage_fu[step] + rho_ullage_fu[step]*V_dot_fu*h_in)
        / (m_ullage_fu[step] + partial_derivative_of_ullage_density_with_respect_to_internal_energy_in_fuel*V_ullage_tank_fu*e_ullage_fu[step] - partial_derivative_of_ullage_density_with_respect_to_internal_energy_in_fuel*V_ullage_tank_fu*h_in)
    ) # [J/kg/s] rate of change of internal energy in fuel ullage
    # Fuel Mdot
    mdot_ullage_fu[step] = partial_derivative_of_ullage_density_with_respect_to_internal_energy_in_fuel*e_dot_fu*V_ullage_tank_fu + rho_ullage_fu[step]*V_dot_fu # [kg/s] rate of change of mass in fuel ullage

    # Advance to next step

    # Update tank properties
    e_ullage_ox[step + 1] = e_ullage_ox[step] + e_dot_ox * DT
    m_ullage_ox[step + 1] = m_ullage_ox[step] + mdot_ullage_ox[step] * DT
    T_ullage_ox[step + 1] = PropsSI('T', 'D', rho_ullage_ox[step], 'U', e_ullage_ox[step + 1], PRESS_GAS)

    e_ullage_fu[step + 1] = e_ullage_fu[step] + e_dot_fu * DT
    m_ullage_fu[step + 1] = m_ullage_fu[step] + mdot_ullage_fu[step] * DT
    T_ullage_fu[step + 1] = PropsSI('T', 'D', rho_ullage_fu[step], 'U', e_ullage_fu[step + 1], PRESS_GAS)

    T_tank_wall_ox[step + 1] = T_tank_wall_ox[step] + (Q_dot_ox_tank * DT) / (m_tank_ox * CP_TANK)
    T_tank_wall_fu[step + 1] = T_tank_wall_fu[step] + (Q_dot_fu_tank * DT) / (m_tank_fu * CP_TANK)

    V_ullage_tank_ox += V_dot_ox * DT
    V_ullage_tank_fu += V_dot_fu * DT

    # Update COPV properties
    mdot_total = mdot_ullage_ox[step] + mdot_ullage_fu[step]
    m_copv[step + 1] = m_copv[step] - mdot_total * DT
    rho_copv[step + 1] = m_copv[step + 1] / V_COPV
    P_copv[step + 1] = PropsSI('P', 'D', rho_copv[step + 1], 'S', s_copv, PRESS_GAS)
    T_copv[step + 1] = PropsSI('T', 'D', rho_copv[step + 1], 'S', s_copv, PRESS_GAS)

    # Time inexorably advances
    step += 1
    time += DT
    #print(f'Ox pressure is {PropsSI('P', 'D', rho_ox, 'U', e_ox[step], PRESS_GAS) * c.PA2PSI:.3f} psi')
    #print(f'Fuel pressure is {PropsSI('P', 'D', rho_fu, 'U', e_fu[step], PRESS_GAS) * c.PA2PSI:.3f} psi')
    #print(f'The COPV pressure is {P_copv[step] * c.PA2PSI:.2f} psi and the tank temperatures are {T_ox[step]:.2f} K in ox and {T_fu[step]:.2f} K in fu at time {time:.3f} s.')
    #print(f'The COPV mass is {rho_copv[step]*V_COPV:.3f} and the tank masses are {rho_ox*V_tank_ox:.3f}  in ox and {rho_fu*V_tank_fu:.3f} in fu for a total of {rho_copv[step]*V_COPV + rho_ox*V_tank_ox + rho_fu*V_tank_fu:.3f}.')

# Set terminal values for plotting
mdot_ullage_ox[0] = mdot_ullage_ox[1]
mdot_ullage_fu[0] = mdot_ullage_fu[1]
mdot_ullage_ox[step] = mdot_ullage_ox[step - 1]
mdot_ullage_fu[step] = mdot_ullage_fu[step - 1]
rho_ullage_ox[step] = rho_ullage_ox[step - 1]
rho_ullage_fu[step] = rho_ullage_fu[step - 1]
h_avg_ullage_ox[step] = h_avg_ullage_ox[step - 1]
h_avg_ullage_fu[step] = h_avg_ullage_fu[step - 1]
h_press_line[0] = h_press_line[1]
Q_dot_press_line[0] = Q_dot_press_line[1]
h_press_line[step] = h_press_line[step - 1]
Q_dot_press_line[step] = Q_dot_press_line[step - 1]

# Get equivalent SCFMs for regulator sizing
scfm_ox = equivalent_SCFM(mdot_ullage_ox, T_copv, PRESS_GAS)
scfm_fu = equivalent_SCFM(mdot_ullage_fu, T_copv, PRESS_GAS)

if TEXT_OUTPUT == True:
    print('RESULTS ------------------------------------------------------------------------------------------------------')
    print(f'The final tank temperatures are {T_ullage_ox[-1]:.3f} [K] in ox and {T_ullage_fu[-1]:.3f} [K] in fuel.')
    print(f'The final COPV properties are {P_copv[-1] * c.PA2PSI:.3f} [psi] and {T_copv[-1]:.3f} [K].')
    print(f'The final masses of gas in the tanks are {m_ullage_ox[-1]:.3f} [kg] in ox and {m_ullage_fu[-1]:.3f} [kg] in fuel.')
    print(f'The average mass flow rate out of COPV during the burn is {(m_copv[0] - m_copv[-1])/drain_time:.3f} [kg/s].')

    # Checks
    # Check ideal gas law
    V_ox = V_ullage_ox + V_dot_ox_nom * drain_time
    Z_ox = P_OX * V_ox / (m_ullage_ox[-1] * R_press_gas * T_ullage_ox[-1])
    V_fu= V_ullage_fu + V_dot_fu_nom * drain_time
    Z_fu = P_FU * V_fu / (m_ullage_fu[-1] * R_press_gas * T_ullage_fu[-1])
    Z_real_ox = PropsSI('Z', 'D', rho_ullage_ox[step], 'U', e_ullage_ox[-1], PRESS_GAS)
    Z_real_fu = PropsSI('Z', 'D', rho_ullage_fu[step], 'U', e_ullage_fu[-1], PRESS_GAS)
    P_real_ox = PropsSI('P', 'D', rho_ullage_ox[step], 'U', e_ullage_ox[-1], PRESS_GAS)
    P_real_fu = PropsSI('P', 'D', rho_ullage_fu[step], 'U', e_ullage_fu[-1], PRESS_GAS)
    print('EQUATION OF STATE --------------------------------------------------------------------------------------------')
    print(f'The real final tank pressures are {P_real_ox * c.PA2PSI:.6f} psi in the oxidizer tank and {P_real_fu * c.PA2PSI:.6f} psi in the fuel tank.')
    print(f'The compressibility factor is {Z_ox:.6f} in the oxidizer tank and {Z_fu:.6f} in the fuel tank.')
    print(f'The compressibility factor should be {Z_real_ox:.6f} in the oxidizer tank and {Z_real_fu:.6f} in the fuel tank.')
    # Check mass balance
    start_mass = PropsSI('D', 'P', P_COPV, 'T', T_COPV, PRESS_GAS) * V_COPV
    end_mass = rho_copv[-1] * V_COPV + m_ullage_ox[-1] + m_ullage_fu[-1]
    print('CONSERVATION OF MASS -----------------------------------------------------------------------------------------')
    print(f'The starting mass of {PRESS_GAS} is {start_mass:.6f} kg and the ending mass is {end_mass:.6f} kg for a percent change of {(end_mass - start_mass)/start_mass * 100:.3f} %')
    # Check energy balance
    start_energy = start_mass * PropsSI('U', 'P', P_COPV, 'T', T_COPV, PRESS_GAS)
    prepress_energy = (
        PropsSI('D', 'P', P0_copv, 'T', T0_copv, PRESS_GAS) * V_COPV * PropsSI('U', 'P', P0_copv, 'T', T0_copv, PRESS_GAS) 
        + (P_OX * V_ullage_ox) / (R_press_gas * T0_ox) * PropsSI('U', 'P', P_OX, 'T', T0_ox, PRESS_GAS)
        + (P_FU * V_ullage_fu) / (R_press_gas * T0_fu) * PropsSI('U', 'P', P_FU, 'T', T0_fu, PRESS_GAS)
    )
    end_energy = rho_copv[-1] * V_COPV * PropsSI('U', 'P', P_copv[-1], 'T', T_copv[-1], PRESS_GAS) + m_ullage_ox[-1] * PropsSI('U', 'P', P_OX, 'T', T_ullage_ox[-1], PRESS_GAS) + m_ullage_fu[-1] * PropsSI('U', 'P', P_FU, 'T', T_ullage_fu[-1], PRESS_GAS)
    work = P_OX * V_dot_ox_nom * drain_time + P_FU * V_dot_fu_nom * drain_time
    net_end_energy = end_energy + work
    print('CONSERVATION OF ENERGY ---------------------------------------------------------------------------------------')
    print(f'The starting energy is {start_energy/1000:.3f} kJ and the ending energy is {end_energy/1000:.3f} kJ.')
    print(f'The PdV work done is {work/1000:.3f} kJ, so the net energy is {net_end_energy/1000:.3f} kJ with a discrepancy of {(net_end_energy - start_energy)/start_energy * 100:.3f} %')
    print(f'The intermediate energy after pre-pressurization is {prepress_energy/1000:.3f} kJ.')

if PLOT_OUTPUT == True:
    
    plt.rcParams['axes.formatter.useoffset'] = False
    plt.rcParams['axes.formatter.use_mathtext'] = False
    plt.rcParams['axes.formatter.limits'] = (-9, 9)
    
    # Plot results
    fig, axs = plt.subplots(2, 3)
    fig.suptitle('Constant-Pressure Pressurization Simulation Results')

    # COPV plots
    axs[0, 0].plot(times, P_copv * c.PA2PSI, 'g')
    axs[0, 0].set_title('COPV Pressure vs. Time')
    axs[0, 0].set_ylabel('Pressure [psi]')
    axs[1, 0].plot(times, T_copv, 'm')
    axs[1, 0].set_title('COPV Temperature vs. Time')
    axs[1, 0].set_ylabel('Temperature [K]')

    # Tank ullage temperatures
    axs[0, 1].plot(times, T_ullage_ox, 'b', times, T_ullage_fu, 'r')
    axs[0, 1].set_title('Tank Ullage Temperature vs. Time')
    axs[0, 1].set_ylabel('Temperature [K]')
    axs[0, 1].legend(['Oxidizer', 'Fuel'])

    if ADIABATIC == False:
        # Heat transfer coefficients
        axs[1, 1].plot(times, h_avg_ullage_ox, 'b', times, h_avg_ullage_fu, 'r')
        axs[1, 1].set_title('Average Ullage Heat Transfer Coefficient vs. Time')
        axs[1, 1].set_ylabel('Heat transfer coefficient [W/m2K]')
        axs[1, 1].legend(['Oxidizer', 'Fuel'])

    # Mass flow rates
    axs[0, 2].plot(times, mdot_ullage_ox, 'b', times, mdot_ullage_fu, 'r')
    axs[0, 2].set_title('Pressurant Mass Flow Rate vs. Time')
    axs[0, 2].set_ylabel('Mass flow rate [kg/s]')
    axs[0, 2].legend(['Oxidizer', 'Fuel'])

    # Equivalent pressurant SCFM
    axs[1, 2].plot(times, scfm_ox, 'b', times, scfm_fu, 'r')
    axs[1, 2].set_title('Pressurant Equivalent N2 SCFM vs. Time')
    axs[1, 2].set_ylabel('Flow rate [SCFM N2]')
    axs[1, 2].legend(['Oxidizer', 'Fuel'])

    for ax in axs.flat:
        ax.set(xlabel='Time [s]')
        ax.grid(True)

    # manager = plt.get_current_fig_manager()
    # manager.window.state('zoomed')
    fig.subplots_adjust(top=0.95, bottom=0.05, hspace=0.3, left=0.05, right=0.95)
    plt.show()


