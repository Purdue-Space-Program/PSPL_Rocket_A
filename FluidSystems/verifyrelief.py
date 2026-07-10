from CoolProp.CoolProp import PropsSI
import numpy as np
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
import matplotlib.pyplot as plt

def mdot_to_SCFM(mdot):
    rho_std = PropsSI('D', 'T', 288.706, 'P', ambient_pressure, 'Nitrogen') # [kg/m^3] 60F standard temp
    qdot_std = mdot / rho_std
    return qdot_std * c.M32FT3 * 60 # [SCFM]

def calc_flow_generant(p_i_psia, dia_in, Kd, temp_rankine):
    # https://www.generant.com/wp-content/uploads/2017/04/Flow_Calculator_Explanation.pdf
    cp = PropsSI('Cpmass', 'T', temp_rankine * c.RANK2KELVIN, 'P', p_i_psia * c.PSI2PA, 'Nitrogen') # [J/kg-K]
    cv = PropsSI('Cvmass', 'T', temp_rankine * c.RANK2KELVIN, 'P', p_i_psia * c.PSI2PA, 'Nitrogen') # [J/kg-K]
    k = cp / cv
    orifice_area_in2 = np.pi * (dia_in / 2)**2
    mdot_lbm_per_second = Kd * orifice_area_in2 * p_i_psia * np.sqrt(((k * 32.174) / (temp_rankine * (1545.35 / 28))) * (2 / (k + 1))**((k + 1) / (k - 1))) # [lbm/s]
    mdot = mdot_lbm_per_second * c.LBM2KG # [kg/s]
    scfm = mdot_to_SCFM(mdot)
    return mdot, scfm

def calc_kd_and_orifice(inlet_size):
    if inlet_size == "1/8":
        dia_in = 0.215
        Kd = 0.57
    elif inlet_size == "1/4" or inlet_size == "3/8":
        dia_in = 0.275
        Kd = 0.65
    elif inlet_size == "1/2":
        dia_in = 0.515
        Kd = 0.35
    else:
        raise ValueError("Invalid inlet size. Choose either 1/8, 1/4, 3/8, or 1/2.")
    return dia_in, Kd

ambient_temp = 293 # [K] Most conservative for flow calcs
ambient_pressure = 1 * c.ATM2PA
set_pressure = 500 # [psia] 
inlet_size = "1/2" # Choose either 1/8, 1/4, 3/8, or 1/2

dia_in, Kd = calc_kd_and_orifice(inlet_size)
mdot, scfm = calc_flow_generant(set_pressure * 1.1, dia_in, Kd, ambient_temp * c.KELVIN2RANK) # Using 110% of nominal set pressure
print("INPUT PARAMETERS")
print(f"Input Pressure (110% of Nominal Set Pressure): {set_pressure * 1.1} psia")
print(f"Input Temperature (Most conservative): {ambient_temp} K")
print(f"Inlet Size: {inlet_size} in")
print("OUTPUT PARAMETERS")
print(f"Calculated Orifice Diameter: {dia_in} in")
print(f"Calculated Discharge Coefficient (Kd): {Kd}")
print(f"Max Mass Flow Rate: {mdot:.2f} kg/s")
print(f"Max SCFM: {scfm:.2f} SCFM")

'''inlet_sizes = ["1/8", "1/4", "3/8", "1/2"]
scfm_values = []
for inlet_size in inlet_sizes:
    dia_in, Kd = calc_kd_and_orifice(inlet_size)
    mdot, scfm = calc_flow_generant(set_pressure * 1.1, dia_in, Kd, ambient_temp * c.KELVIN2RANK)
    scfm_values.append(scfm)

plt.plot(inlet_sizes, scfm_values, marker='o')
plt.xlabel('Inlet Size (in)')
plt.ylabel('Max SCFM')
plt.title('Flow Rate vs Inlet Size')
plt.show()'''

pressure_values = np.linspace(400, 1000, 100) # [psia]
scfm_pressure_values = []
for p in pressure_values:
    dia_in, Kd = calc_kd_and_orifice("1/2")
    mdot, scfm = calc_flow_generant(p * 1.1, dia_in, Kd, ambient_temp * c.KELVIN2RANK)
    scfm_pressure_values.append(scfm)

plt.plot(pressure_values, scfm_pressure_values)
plt.xlabel('Pressure (psia)')
plt.ylabel('Max SCFM')
plt.title('Flow Rate vs Pressure')
plt.show()