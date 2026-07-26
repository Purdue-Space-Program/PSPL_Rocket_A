from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt
import sys,pathlib,collections,importlib.abc; r=next(p for p in pathlib.Path(__file__).resolve().parents if p.name=="PSPL_Rocket_A");   m=collections.defaultdict(list); [m[f.stem.casefold()].append(f) for f in r.rglob("*.py") if f.stem.isidentifier() and f.name != "__init__.py" and not (set(f.relative_to(r).parts)&{".git",".venv","__pycache__","build","dist"})]; dup={k:v for k,v in m.items() if   len(v)>1}; sys.meta_path.insert(0,type("AmbiguousBareImportBlocker",(importlib.abc.MetaPathFinder,),{"find_spec":lambda   self,fullname,path=None,target=None: (_ for _ in ()).throw(ImportError(f"The import {fullname!r} could refer to any following packages: "+" ".join("\n\t" + str(p.relative_to(r)) for p in dup[fullname.casefold()])+f"\n\nSpecify which package it is by using the folder.\nFor example:\n\t'import SFD.{fullname}'")) if "." not in fullname and fullname.casefold() in dup else None})()); sys.path.insert(0,str(r)) if str(r) not   in sys.path else None; [sys.path.append(str(v[0].parent)) for k,v in m.items() if k not in dup and str(v[0].parent) not in sys.path]
import constants as c

def mdot_to_SCFM(mdot):
    rho_std = PropsSI('D', 'T', 288.706, 'P', ambient_pressure, 'Nitrogen') # [kg/m^3] 60F standard temp
    qdot_std = mdot / rho_std
    return qdot_std * c.M32FT3 * 60 # [SCFM]

def calc_flow_generant(inlet_pressure_psia, dia_in, Kd, temp_rankine):
    # https://www.generant.com/wp-content/uploads/2017/04/Flow_Calculator_Explanation.pdf
    cp = PropsSI('Cpmass', 'T', temp_rankine * c.RANK2KELVIN, 'P', inlet_pressure_psia * c.PSI2PA, 'Nitrogen') # [J/kg-K]
    cv = PropsSI('Cvmass', 'T', temp_rankine * c.RANK2KELVIN, 'P', inlet_pressure_psia * c.PSI2PA, 'Nitrogen') # [J/kg-K]
    k = cp / cv
    orifice_area_in2 = np.pi * (dia_in / 2)**2
    mdot_lbm_per_second = Kd * orifice_area_in2 * inlet_pressure_psia * np.sqrt(((k * 32.174) / (temp_rankine * (1545.35 / 28))) * (2 / (k + 1))**((k + 1) / (k - 1))) # [lbm/s]
    mdot = mdot_lbm_per_second * c.LBM2KG # [kg/s]
    scfm = mdot_to_SCFM(mdot)
    return mdot, scfm

def calc_kd_and_orifice(inlet_size): # these values are for HRPVA: https://www.generant.com/wp-content/uploads/2024/12/SA.SL_.HPRV001_H.4864-HPRV-Series-Product-Literature.pdf
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

ambient_temp = 300 # [K]
ambient_pressure = 1 * c.ATM2PA
relief_set_pressure = 550 + (1 * c.ATM2PA * c.PA2PSI) # [psia]: https://www.globaltestsupply.com/product/generant-hprva-series-vent-to-atmosphere-high-pressure-relief-valve?cfrconfigurator_id=F172453&cfrkey=HPRVA-500B-T-500&cfrmodel=500B&cfrseal=T&cfrcrack=500 
relief_valve_overpressure_factor = 1.1 # 10% overpressure for flow calcs
relief_valve_inlet_size = "1/2" # Choose either 1/8, 1/4, 3/8, or 1/2

relief_valve_inlet_diameter_inches, Kd = calc_kd_and_orifice(relief_valve_inlet_size)
relief_valve_mass_flow_rate, relief_valve_scfm = calc_flow_generant(relief_set_pressure * relief_valve_overpressure_factor, relief_valve_inlet_diameter_inches, Kd, ambient_temp * c.KELVIN2RANK) # Using 110% of nominal set pressure
print("INPUT PARAMETERS")
print(f"\tInlet Pressure (110% of Nominal Set Pressure): {relief_set_pressure * relief_valve_overpressure_factor:.2f} psia")
print(f"\tInlet Temperature (Most conservative): {ambient_temp} K")
print(f"\tInlet Size: {relief_valve_inlet_size} in")
print("OUTPUT PARAMETERS")
print(f"\tCalculated Orifice Diameter: {relief_valve_inlet_diameter_inches} in")
print(f"\tCalculated Discharge Coefficient (Kd): {Kd}")
print(f"\tMax Mass Flow Rate: {relief_valve_mass_flow_rate:.2f} kg/s")
print(f"\tMax SCFM: {relief_valve_scfm:.2f} SCFM")

print_all_relief_sizes = False
if print_all_relief_sizes:
    inlet_sizes = ["1/8", "1/4", "3/8", "1/2"]
    scfm_values = []
    for inlet_size in inlet_sizes:
        dia_in, Kd = calc_kd_and_orifice(inlet_size)
        relief_valve_mass_flow_rate, relief_valve_scfm = calc_flow_generant(relief_set_pressure * relief_valve_overpressure_factor, dia_in, Kd, ambient_temp * c.KELVIN2RANK)
        scfm_values.append(relief_valve_scfm)

    plt.plot(inlet_sizes, scfm_values, 'o')
    plt.xlabel('Inlet Size (in)')
    plt.ylabel('Max SCFM')
    plt.title('Flow Rate vs Inlet Size')
    plt.show()

pressure_values = np.arange(400 + (1 * c.ATM2PA * c.PA2PSI), 1000 + (1 * c.ATM2PA * c.PA2PSI), 50) # sold in increments of 50 PSIG [PSIA]
scfm_pressure_values = []
for pressure in pressure_values:
    relief_valve_inlet_diameter_inches, Kd = calc_kd_and_orifice("1/2")
    relief_valve_mass_flow_rate, relief_valve_scfm = calc_flow_generant(pressure * relief_valve_overpressure_factor, relief_valve_inlet_diameter_inches, Kd, ambient_temp * c.KELVIN2RANK)
    scfm_pressure_values.append(relief_valve_scfm)

plt.plot(pressure_values, scfm_pressure_values, 'o')
plt.xlabel('Relief Valve Inlet Pressure [PSIA]')
plt.ylabel('Max SCFM')
plt.title('Flow Rate vs Pressure')
plt.show()