# Constant-Pressure Press Sim Functions
# These functions go with the constant-pressure press sim
# Hugo Filmer

from CoolProp.CoolProp import PropsSI
import numpy as np

def adiabatic_press(P1, T1, V_COPV, P_tank1, V_tank1, P_tank2, V_tank2, gas):
    # Models simultaneous adiabatic pressurization of two tanks from an isentropic COPV

    # Constants
    STOP_DIFF = 0.001

    # COPV starting conditions
    s1 = PropsSI('S', 'P', P1, 'T', T1, gas) # [J/kgK] The starting COPV entropy
    rho1 = PropsSI('D', 'P', P1, 'T', T1, gas) # [kg/m^3] The starting COPV density
    e1 = PropsSI('U', 'P', P1, 'T', T1, gas) # [J/kg] The starting COPV specific internal energy
    m1 = rho1 * V_COPV # [kg] The starting COPV mass
    E1 = m1 * e1 # [J] The starting COPV energy

    # Iterative process
    T_guess = 400 # [K] Guessed final tank temperature to start
    T_tank1 = T_guess # [K] The final tank temperature for tank 1
    T_tank2 = T_guess # [K] The final tank temperature for tank 2
    T_iter1 = 0 # [K] The final tank temperature for tank 1 (variable for iteration)
    T_iter2 = 0 # [K] The final tank temperature for tank 2 (variable for iteration)
    iter_count = 0 # [] The number of iterations

    while abs((T_iter1 - T_tank1)/T_tank1)*100 > STOP_DIFF and abs((T_iter2 - T_tank2)/T_tank2)*100 > STOP_DIFF: # While percent difference greater than

        T_iter1 = T_tank1
        T_iter2 = T_tank2

        iter_count += 1

        # The final tank mass is the sum of tank density * tank volume
        rho_tank1 = PropsSI('D', 'P', P_tank1, 'T', T_iter1, gas)
        rho_tank2 = PropsSI('D', 'P', P_tank2, 'T', T_iter2, gas)
        m_tank1 = V_tank1 * rho_tank1
        m_tank2 = V_tank2 * rho_tank2

        # Mass is conserved, so the final COPV mass is the starting COPV mass - the mass in the tanks
        m2_copv = m1 - m_tank1 - m_tank2

        # Density is defined as mass / volume
        rho2_copv = m2_copv / V_COPV

        # The expansion of gas in the COPV is isentropic
        s2 = s1

        # The final COPV internal energy is specific internal energy * mass
        e2_copv = PropsSI('U', 'D', rho2_copv, 'S', s2, gas)
        E2_copv = m2_copv * e2_copv

        # Energy is conserved, so the final tank energy is the starting COPV energy - the final COPV energy
        E_tanks = E1 - E2_copv
        
        # Assume energy is divided with pressure * volume (ideal gas model)
        E_tank1 = E_tanks * (P_tank1 * V_tank1) / (P_tank1 * V_tank1 + P_tank2 * V_tank2)
        E_tank2 = E_tanks * (P_tank2 * V_tank2) / (P_tank1 * V_tank1 + P_tank2 * V_tank2)

        # Specific internal energy is defined as energy / mass
        e_tank1 = E_tank1 / m_tank1
        e_tank2 = E_tank2 / m_tank2

        # Calculate the new final tank temperature estimate
        T_tank1 = PropsSI('T', 'P', P_tank1, 'U', e_tank1, gas)
        T_tank2 = PropsSI('T', 'P', P_tank2, 'U', e_tank2, gas)

    # Calculate other final COPV parameters
    T2 = PropsSI('T', 'D', rho2_copv, 'S', s2, gas)
    P2 = PropsSI('P', 'D', rho2_copv, 'S', s2, gas)

    return [P2, T2, T_tank1, T_tank2, iter_count]

def isothermal_press(P1, T1, V_COPV, P_tank1, V_tank1, T_tank1, P_tank2, V_tank2, T_tank2, gas):
    # Models prepressurization of two tanks from an isentropic COPV with known tank temperatures

    # COPV starting conditions
    s1 = PropsSI('S', 'P', P1, 'T', T1, gas) # [J/kgK] The starting COPV entropy
    rho1 = PropsSI('D', 'P', P1, 'T', T1, gas) # [kg/m^3] The starting COPV density
    m1 = rho1 * V_COPV # [kg] The starting COPV mass

     # The final tank mass is the sum of tank density * tank volume
    rho_tank1 = PropsSI('D', 'P', P_tank1, 'T', T_tank1, gas)
    rho_tank2 = PropsSI('D', 'P', P_tank2, 'T', T_tank2, gas)
    m_tank1 = V_tank1 * rho_tank1
    m_tank2 = V_tank2 * rho_tank2

    # Mass is conserved, so the final COPV mass is the starting COPV mass - the mass in the tanks
    m2_copv = m1 - m_tank1 - m_tank2

    # Density is defined as mass / volume
    rho2_copv = m2_copv / V_COPV

    # The expansion of gas in the COPV is isentropic
    s2 = s1

    # Calculate other final COPV parameters
    T2 = PropsSI('T', 'D', rho2_copv, 'S', s2, gas)
    P2 = PropsSI('P', 'D', rho2_copv, 'S', s2, gas)

    return [P2, T2, T_tank1 + 0.001, T_tank2 + 0.001, 0]

def wall_convect(g, P_tank, T_s, T_inf, L, fluid):

    T_film = (T_s + T_inf)/2

    rho = PropsSI('D', 'P', P_tank, 'T', T_film, fluid)
    nu = PropsSI('V', 'P', P_tank, 'T', T_film, fluid) / rho
    beta = -1/rho * PropsSI('d(D)/d(T)|P', 'P', P_tank, 'T', T_film, fluid)

    Gr = abs((g * beta * (T_s - T_inf) * L**3) / nu**2)
    Pr = PropsSI('PRANDTL', 'P', P_tank, 'T', T_film, fluid) / PropsSI('D', 'P', P_tank, 'T', T_film, fluid)
    Ra = Gr * Pr

    Nu = 0.1 * Ra**(1/3)
    k = PropsSI('L', 'P', P_tank, 'T', T_film, fluid)
    h = Nu * k / L

    return h

def plate_convect(g, P_tank, T_s, T_inf, L, fluid, enhancement, film=True):

    if film == True:
        T_film = (T_s + T_inf)/2
    else:
        T_film = T_inf

    rho = PropsSI('D', 'P', P_tank, 'T', T_film, fluid)
    nu = PropsSI('V', 'P', P_tank, 'T', T_film, fluid) / rho
    beta = -1/rho * PropsSI('d(D)/d(T)|P', 'P', P_tank, 'T', T_film, fluid)

    Gr = abs((g * beta * (T_s - T_inf) * L**3) / nu**2)
    Pr = PropsSI('PRANDTL', 'P', P_tank, 'T', T_film, fluid) / PropsSI('D', 'P', P_tank, 'T', T_film, fluid)
    Ra = Gr * Pr

    if enhancement == 'enhanced':
        Nu = 0.1 * Ra**(1/3)
    elif enhancement == 'reduced':
        Nu = 0.27 * Ra**(1/4)

    k = PropsSI('L', 'P', P_tank, 'T', T_film, fluid)
    h = Nu * k / L

    return h

def interface_convect(g, P_tank, T1, T2, L, fluid1, fluid2):

    T_interface = (T1 + T2)/2 # Initial guess

    for i in range(5):

        # Upper surface
        if T1 > T2:
            enhancement = 'reduced'
        else:
            enhancement = 'enhanced'
        h_upper = plate_convect(g, P_tank, T_interface, T1, L, fluid1, enhancement)

        # Lower surface
        if T1 > T2:
            enhancement = 'enhanced'
        else:
            enhancement = 'reduced'
        h_lower = plate_convect(g, P_tank, T_interface, T2, L, fluid2, enhancement, film=False)

        T_interface = (h_upper * T1 + h_lower * T2) / (h_upper + h_lower)

    q_interface = h_upper * (T1 - T_interface)
    h_effective = q_interface / (T1 - T2)

    return h_effective

def line_heat_transfer(P_in, h_in, T_line, D_line, L_line, m_dot, gas):

    # Dimensions
    e_line = 0.0015 / 1000 # [m] line relative roughness (typical for drawn tubing)
    A_line = np.pi * D_line / 4
    As_line = np.pi * L_line * D_line

    # Evaluate properties
    T_in = PropsSI('T', 'P', P_in, 'H', h_in, gas)
    k = PropsSI('L', 'P', P_in, 'H', h_in, gas)
    Cp = PropsSI('C', 'P', P_in, 'H', h_in, gas)

    # Estimate heat transfer coefficient
    Nu = 3.66
    h = Nu * k / D_line

    # Calculate temperatures
    T_out = T_line - (T_line - T_in)*np.exp(-h*As_line / (m_dot * Cp))
    delta_T_lm = (T_in - T_out) / np.log((T_line - T_out) / (T_line - T_in))

    # Calculate heat transfer
    Q_dot = h * As_line * delta_T_lm

    return [h, Q_dot]

def dZ_dT(P, T, gas):
    Z = PropsSI('Z', 'P', P, 'T', T, gas)
    rho = PropsSI('D', 'P', P, 'T', T, gas)
    drho_dT = PropsSI('d(D)/d(T)|P', 'P', P, 'T', T, gas)

    dZ_dT = -Z*drho_dT/rho - Z/T

    return dZ_dT

def equivalent_SCFM(mdots, Ts, gas):
    # Determines equivalent N2 SCFM for regulator sizing

    # Conversion factors
    M32FT3 = 35.3147
    SEC2MIN = 60

    if gas == 'helium':
        Fg = 2.65 # [] gas correction factor
    elif gas == 'nitrogen':
        Fg = 1

    Ft = -0.00212133*Ts + 1.61586 # [] temperature correction factor based on linear regression from Swagelok data

    rho_ambient = PropsSI('D', 'P', 101325, 'T', 293, gas)
    scfm = mdots / rho_ambient * M32FT3 * SEC2MIN

    corrected_scfm = scfm / (Fg * Ft)

    return corrected_scfm