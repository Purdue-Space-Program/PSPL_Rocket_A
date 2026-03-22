from CoolProp.CoolProp import PropsSI
import numpy as np
from pint import UnitRegistry
import constants as c

u = UnitRegistry()

nitrogen_standard_conditions_density = PropsSI("D", "P", 1 * c.ATM2PA, "T", c.T_AMBIENT, "nitrogen") * (u.kilogram / (u.meter**3))
print(f"\nnitrogen_standard_conditions_density: {nitrogen_standard_conditions_density:.4g}\n")

pressure_vessel_fluid_pressure = 500 * c.PSI2PA * u.pascal
pressure_vessel_fluid_temperature = c.T_AMBIENT * u.kelvin
pressure_vessel_fluid_density = PropsSI("D", "P", pressure_vessel_fluid_pressure.magnitude, "T", pressure_vessel_fluid_temperature.magnitude, "nitrogen") * (u.kilogram / (u.meter**3))
pressure_vessel_orifice_area = np.pi * (((0.050 * c.IN2M) / 2)**2) * (u.meter**2)
pressure_vessel_fluid_pressure_drop = pressure_vessel_fluid_pressure - 1 * c.ATM2PA * u.pascal

def CalculateMassFlowRate(area, Cd, rho, pressure_drop):
    m_dot = area * (Cd * np.sqrt(2 * rho * pressure_drop))
    m_dot = m_dot.to(u.kilogram / u.second)
    return m_dot


pressure_vessel_fluid_mass_flow_rate = CalculateMassFlowRate(pressure_vessel_orifice_area,
                                                             0.6,
                                                             pressure_vessel_fluid_density,
                                                             pressure_vessel_fluid_pressure_drop,
                                                             )


print(f"pressure_vessel_fluid_mass_flow_rate: {pressure_vessel_fluid_mass_flow_rate:.4g} kg/s")
print(f"pressure_vessel_fluid_mass_flow_rate: {pressure_vessel_fluid_mass_flow_rate * c.KG2G:.4g} g/s\n")

SCFM_2_KG_S = ((1**3 * c.FT32M3) * nitrogen_standard_conditions_density) / 60

print(f"\nSCFM2KG_S: {SCFM_2_KG_S:.4g} kg/s")
print(f"SCFM2KG_S: {SCFM_2_KG_S * c.KG2G:.4g} g/s\n")
