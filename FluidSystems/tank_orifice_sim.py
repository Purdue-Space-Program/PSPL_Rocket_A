import numpy as np
import sys
import os
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt
from dataclasses import dataclass
import copy

os.chdir(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import constants as c
import Press_Cv_Calcs
import vehicle_parameters

parameters, wet_mass_distribution, dry_mass_distribution = vehicle_parameters.main()

def Cv_Choked_from_SCFM(press_gas, SCFM_pressurant, P1, T1):

    # Unit conversions
    SEC2MIN = 1 / 60
    K2R = 9 / 5

    # Constants
    N2 = 22.67 # units constant for flow rate equation

    # Standard conditions
    P_STD = c.ATM2PA # [pa] standard pressure
    T_STD = 288.7 # [K] standard temperature

    rho_std_air = PropsSI('D', 'P', P_STD, 'T', T_STD, 'air') # [kg/^3] standard density of air
    rho_std_press_gas = PropsSI('D', 'P', P_STD, 'T', T_STD, press_gas) # [kg/m^3] standard density of press gas

    Gg_pressurant = rho_std_press_gas / rho_std_air # [] specific gravity of pressurant

    P1_psi = P1 * c.PA2PSI # [psia] inlet pressure in psia
    T1_r = T1 * K2R # [deg R] inlet temperature in deg R

    Cv = SCFM_pressurant * np.sqrt(Gg_pressurant * T1_r) / (0.471 * N2 * P1_psi)

    return Cv

def Calculate_Choked_Mass_Flow_from_Cv(press_gas, Cv, P1, T1):

    # Unit conversions
    SEC2MIN = 1 / 60
    K2R = 9 / 5

    # Constants
    N2 = 22.67 # units constant for flow rate equation

    # Standard conditions
    P_STD = c.ATM2PA # [pa] standard pressure
    T_STD = 288.7 # [K] standard temperature


    rho_std_air = PropsSI('D', 'P', P_STD, 'T', T_STD, 'air') # [kg/^3] standard density of air
    rho_std_press_gas = PropsSI('D', 'P', P_STD, 'T', T_STD, press_gas) # [kg/m^3] standard density of press gas

    Gg_press_gas = rho_std_press_gas / rho_std_air # [] specific gravity of press gas

    P1_psi = P1 * c.PA2PSI # [psia] inlet pressure in psia
    T1_r = T1 * K2R # [deg R] inlet temperature in deg R

    SCFM_press_gas = Cv / (np.sqrt(Gg_press_gas * T1_r) / (0.471 * N2 * P1_psi))

    q_std_press_gas = SCFM_press_gas / (c.M32FT3 / SEC2MIN) # [m^3/s] standard volumetric flow rate of press gas

    m_dot = q_std_press_gas * rho_std_press_gas

    return m_dot


def Simulate_Orifice_Emptying(tank_with_orifice):

    fig, ax = plt.subplots(2, 2, sharex=False)

    tank_with_orifice.orifice.flow_coefficient = Cv_Choked_from_SCFM(tank_with_orifice.pressurant_name, tank_with_orifice.orifice.test_flow_rate, tank_with_orifice.orifice.test_pressure, 288.7)
    print(f"orifice {tank_with_orifice.orifice.nominal_size}: {tank_with_orifice.orifice.flow_coefficient:.5f} Cv")


    dt = 0.5

    for model_type in tank_with_orifice.model_types:

        initial_time = 0
        pressurant_initial_mass_flow_rate = 0 # [kg/s]
        pressurant_initial_temperature = tank_with_orifice.pressurant_initial_temperature # [k]
        pressurant_initial_pressure = tank_with_orifice.pressurant_initial_pressure
        if tank_with_orifice.pressurant_name == "oxygen":
            pressurant_initial_density = PropsSI("D", "P", pressurant_initial_pressure, "Q", 1, tank_with_orifice.pressurant_name)
            pressurant_initial_entropy = PropsSI("S", "P", pressurant_initial_pressure, "Q", 1, tank_with_orifice.pressurant_name)
            final_pressure = 30
        elif tank_with_orifice.pressurant_name == "nitrogen":
            pressurant_initial_density = PropsSI("D", "T", pressurant_initial_temperature, "P", pressurant_initial_pressure, tank_with_orifice.pressurant_name)
            pressurant_initial_entropy = PropsSI("S", "T", pressurant_initial_temperature, "P", pressurant_initial_pressure, tank_with_orifice.pressurant_name)
            final_pressure = 100

        pressurant_initial_mass = tank_with_orifice.volume * pressurant_initial_density

        current_time = initial_time
        pressurant_current_mass_flow_rate = pressurant_initial_mass_flow_rate
        pressurant_current_density = pressurant_initial_density
        pressurant_current_pressure = pressurant_initial_pressure
        pressurant_current_mass = pressurant_initial_mass
        pressurant_current_temperature = pressurant_initial_temperature


        time_array = np.array([], dtype=float)
        pressurant_mass_flow_rate_array = np.array([], dtype=float)
        pressurant_pressure_array = np.array([], dtype=float)
        pressurant_density_array = np.array([], dtype=float)
        pressurant_mass_array = np.array([], dtype=float)
        pressurant_temperature_array = np.array([], dtype=float)

        time_array = np.append(time_array, current_time)
        pressurant_mass_flow_rate_array = np.append(pressurant_mass_flow_rate_array, pressurant_current_mass_flow_rate)
        pressurant_pressure_array = np.append(pressurant_pressure_array, pressurant_current_pressure)
        pressurant_density_array = np.append(pressurant_density_array, pressurant_current_density)
        pressurant_mass_array = np.append(pressurant_mass_array, pressurant_current_mass)
        pressurant_temperature_array = np.append(pressurant_temperature_array, pressurant_current_temperature)


        while pressurant_current_pressure > (final_pressure * c.PSI2PA):
            # print(f"pressurant_current_pressure: {pressurant_current_pressure}")
            pressurant_current_mass_flow_rate = Calculate_Choked_Mass_Flow_from_Cv(tank_with_orifice.pressurant_name, tank_with_orifice.orifice.flow_coefficient, pressurant_current_pressure, pressurant_current_temperature)

            new_time = time_array[-1] + dt
            pressurant_new_mass = pressurant_mass_array[-1] - (pressurant_current_mass_flow_rate * dt)
            pressurant_new_density = pressurant_new_mass / tank_with_orifice.volume

            if model_type == "isentropic expansion":
                pressurant_new_entropy = pressurant_initial_entropy
                pressurant_new_temperature = PropsSI("T", "S", pressurant_new_entropy, "D", pressurant_new_density, tank_with_orifice.pressurant_name)
                pressurant_new_pressure = PropsSI("P", "S", pressurant_new_entropy, "D", pressurant_new_density, tank_with_orifice.pressurant_name)
            elif model_type == "isothermal expansion":
                pressurant_new_temperature = pressurant_current_temperature
                pressurant_new_pressure = PropsSI("P", "T", pressurant_new_temperature, "D", pressurant_new_density, tank_with_orifice.pressurant_name)
            elif model_type == "boil-off":
                pressurant_new_temperature = pressurant_current_temperature
                pressurant_new_pressure = pressurant_current_pressure

            quality_value = PropsSI("Q", "T", pressurant_new_temperature, "D", pressurant_new_density, tank_with_orifice.pressurant_name)
            # print("Q =", quality_value)

            time_array = np.append(time_array, new_time)
            pressurant_mass_flow_rate_array = np.append(pressurant_mass_flow_rate_array, pressurant_current_mass_flow_rate)
            pressurant_pressure_array = np.append(pressurant_pressure_array, pressurant_new_pressure)
            pressurant_density_array = np.append(pressurant_density_array, pressurant_new_density)
            pressurant_mass_array = np.append(pressurant_mass_array, pressurant_new_mass)
            pressurant_temperature_array = np.append(pressurant_temperature_array, pressurant_new_temperature)

            current_time = new_time
            pressurant_current_mass = pressurant_new_mass
            pressurant_current_density = pressurant_new_density
            # pressurant_current_entropy = pressurant_new_entropy
            pressurant_current_temperature = pressurant_new_temperature
            pressurant_current_pressure = pressurant_new_pressure

            # print(f"pressurant_current_pressure: {pressurant_current_pressure * c.PA2PSI:.2f}")



        if max(time_array) < 120:
            time_unit_name = "Seconds"
        else:
            time_array = time_array/60
            time_unit_name = "Minutes"

        ax[0,0].plot(time_array, pressurant_pressure_array * c.PA2PSI)
        ax[0,0].set_ylabel(f"Pressure [PSI]")
        ax[0,0].set_xlabel(f"Time [{time_unit_name}]")
        ax[0,0].grid()

        ax[0,1].plot(time_array, pressurant_mass_flow_rate_array)
        ax[0,1].set_ylabel(f"Mass flow rate [kg/s]")
        ax[0,1].set_xlabel(f"Time [{time_unit_name}]")
        ax[0,1].grid()

        ax[1,0].plot(time_array, pressurant_temperature_array)
        ax[1,0].set_ylabel(f"Temperature [K]")
        ax[1,0].set_xlabel(f"Time [{time_unit_name}]")
        ax[1,0].grid()

        ax[1,1].plot(time_array, pressurant_density_array)
        ax[1,1].set_ylabel(f"Density [kg/m^3]")
        ax[1,1].set_xlabel(f"Time [{time_unit_name}]")
        ax[1,1].grid()

        fig.suptitle(f"Emptying {tank_with_orifice.pressurant_name} in {tank_with_orifice.name} through {tank_with_orifice.orifice.nominal_size} orifice using {model_type} model")
        # plt.xlabel("inlet pressure [psi]")
        # plt.ylabel("mass flow rate [kg]")

        fig.legend(tank_with_orifice.model_types)
    plt.show()

def main():

    @dataclass
    class Orifice:
        nominal_size: str
        test_flow_rate: float
        test_pressure: float
        flow_coefficient: float = None

    @dataclass
    class TankWithOrifice:
        name: str
        pressurant_initial_pressure: float
        pressurant_initial_temperature: float
        volume: float
        orifice: Orifice
        model_types: str
        pressurant_name: str

    # data from: https://catalog.okeefecontrols.com/Download/Standard-Conditions-in-Flow-Measurement-5-18-21.pdf
    orifice_0010 = Orifice(
                            nominal_size = "0.010",
                            test_flow_rate = 8.12/60, # value from website is in SCFH [SCFM]
                            test_pressure = (80 * c.PSI2PA) + c.ATM2PA,
                          )

    # data from: https://catalog.okeefecontrols.com/Download/Standard-Conditions-in-Flow-Measurement-5-18-21.pdf
    orifice_0125 = Orifice(
                            nominal_size = "0.125",
                            test_flow_rate = 1263/60, # value from website is in SCFH [SCFM]
                            test_pressure = (80 * c.PSI2PA) + c.ATM2PA,
                          )

    # orifice_0.... = Orifice(
    #                         nominal_size = "0.125",
    #                         test_flow_rate = 1263/60, # [SCFM]
    #                         test_pressure = (80 * c.PSI2PA) + c.ATM2PA,
    #                       )

    COPV_with_orifice = TankWithOrifice(
                                        name = "COPV",
                                        pressurant_initial_pressure = parameters.COPV_starting_pressure,
                                        pressurant_initial_temperature = c.T_AMBIENT,
                                        volume = parameters.COPV_volume,
                                        orifice = orifice_0010,
                                        model_types = ["isothermal expansion", "isentropic expansion"],
                                        pressurant_name = parameters.pressurant_name,
                                       )

    oxidizer_tank_with_orifice = TankWithOrifice(
                                        name = "Liquid Oxygen Tank",
                                        pressurant_initial_pressure = 250 * c.PSI2PA,
                                        pressurant_initial_temperature = PropsSI("T", "P", 50 * c.PSI2PA, "Q", 1, "Oxygen"),
                                        volume = parameters.oxidizer_tank_usable_volume,
                                        orifice = orifice_0010,
                                        model_types = ["isothermal expansion"],
                                        pressurant_name = parameters.oxidizer_name,
                                       )

    Simulate_Orifice_Emptying(COPV_with_orifice)
    # Simulate_Orifice_Emptying(oxidizer_tank_with_orifice)


if __name__ == "__main__":
    main()
