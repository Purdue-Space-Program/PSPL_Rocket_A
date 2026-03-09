from CoolProp.CoolProp import PropsSI
import numpy as np
import final_pathfinder_press_sim
import constants as c
import press_sim

def CvCrit(qdot, Gs, T, P1):
    '''
    qdot: mass flow rate (ft^3/min)
    Gs: Specific gravity of gas (relative to specific_gravity_of_air=1)
    T: Absolute Upstream Temperature (Rankine)
    P1: Inlet pressure (psia)
    cv: Flow coefficient
    '''
    cv = qdot / (0.471 * 22.67 * P1 * np.sqrt(1/(Gs * T))) # https://www.swagelok.com/downloads/webcatalogs/en/ms-06-84.pdf (page 3)
    return cv

# COPV
copv_temp = final_pathfinder_press_sim.T_COPV * c.KELVIN2RANK # [R]
copv_pressure = 4500 # final_pathfinder_press_sim.P_COPV * c.PA2PSI # [psia]
copv_density = PropsSI('D', 'T', copv_temp / c.KELVIN2RANK, 'P', copv_pressure / c.PA2PSI, 'Nitrogen') # kg/m^3
air_density = PropsSI('D', 'T', 293.15, 'P', 101325, 'Air') # kg/m^3
copv_gs = copv_density / air_density # Specific gravity of gas (relative to specific_gravity_of_air=1)
# copv_qdot = final_pathfinder_press_sim.maximum_regulator_mass_flow_rate / copv_density * c.M32FT3 * 60 # [ft^3/min]
regulator_cv = 0.8 # CvCrit(copv_qdot, copv_gs, copv_temp, copv_pressure)
qdot_calculated = 0.471 * 22.67 * regulator_cv * copv_pressure * np.sqrt(1/(copv_gs * copv_temp)) # [ft^3/min]
scfm_copv = qdot_calculated * copv_pressure / 14.7 * 520 / copv_temp # [scfm]

print(f"COPV Temperature: {copv_temp:.2f} R")
print(f"COPV Pressure: {copv_pressure:.2f} psia")
print(f"COPV Density: {copv_density:.2f} kg/m^3")
print(f"Air Density at COPV Conditions: {air_density:.2f} kg/m^3")
print(f"COPV Specific Gravity: {copv_gs:.2f}")
# print(f"COPV Volumetric Flow Rate: {copv_qdot:.2f} ft^3/min")
print(f"COPV Cv: {regulator_cv:.10f}")
print(f"Calculated Volumetric Flow Rate: {qdot_calculated:.2f} ft^3/min")
print(f"Calculated SCFM: {scfm_copv:.2f} scfm")

# LOx Relief Valve
ox_pressure = final_pathfinder_press_sim.P_OX * c.PA2PSI # [psia]
ox_temp = 293 * c.KELVIN2RANK # final_pathfinder_press_sim.T_fill_ox * c.KELVIN2RANK # [R]
ox_tank_density_pressurant = PropsSI('D', 'T', ox_temp / c.KELVIN2RANK, 'P', ox_pressure / c.PA2PSI, 'Nitrogen') # kg/m^3
ox_tank_density_air = PropsSI('D', 'T', 293.15, 'P', 101325, 'Air') # kg/m^3
ox_tank_gs = ox_tank_density_pressurant / ox_tank_density_air # Specific gravity of gas (relative to specific_gravity_of_air=1)
ox_relief_valve_cv = CvCrit(qdot_calculated, ox_tank_gs, ox_temp, ox_pressure)
print(f"LOx Pressure: {ox_pressure:.2f} psia")
print(f"LOx Temperature: {ox_temp:.2f} R")
print(f"LOx Tank Density (Pressurant): {ox_tank_density_pressurant:.2f} kg/m^3")
print(f"LOx Tank Density (Air): {ox_tank_density_air:.2f} kg/m^3")
print(f"LOx Tank Specific Gravity: {ox_tank_gs:.2f}")
print(f"LOx Relief Valve Cv: {ox_relief_valve_cv:.10f}")

'''
# Fuel Relief Valve
fuel_pressure = final_pathfinder_press_sim.P_FU * c.PA2PSI # [psia]
fuel_temp = 293 * c.KELVIN2RANK # final_pathfinder_press_sim.T_fill_ox * c.KELVIN2RANK # [R]
fuel_relief_valve_cv = CvCrit(qdot_calculated, copv_gs, fuel_temp, fuel_pressure)
print(f"Fuel Relief Valve Cv: {fuel_relief_valve_cv:.10f}")
'''