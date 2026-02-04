# This code reflects the vehicle parameters as stated on the vehicle parameters page: https://purdue-space-program.atlassian.net/wiki/x/Aw5qYw
# Import this file into your script as a way to easily and reliably pull parameters. 
# Author: David Gustafsson

from dataclasses import dataclass, fields, field
from typing import Iterator
import matplotlib.pyplot as plt
import numpy as np
import csv
from pathlib import Path
import inspect
import subprocess
from pathlib import Path
import sys
from datetime import datetime
import io

import constants as c
import vehicle_parameters_functions

@dataclass
class VehicleParameters:

    # Structural Parameters
    yield_FoS: float = 1.5
    ultimate_FoS: float = 2.0
    proof_factor: float = 1.5
    bending_FoS: float = None

    # General Parameters
    fuel_name: str = "isopropyl alcohol" 
    oxidizer_name: str = "liquid oxygen"
    tube_outer_diameter: float = 6.0 * c.IN2M         # Outer diameter of tube used in some sections of the rocket
    tube_inner_diameter: float = 5.75 * c.IN2M        # Inner diameter of tube used in some sections of the rocket
    
    # Engine Parameters
    chamber_pressure: float = 250 * c.PSI2PA                                    # The target combustion pressure in the engine [Pascals]
    jet_thrust: float = 668.0 * c.LBF2N                                         # The targeted engine thrust (not accounting for exhaust gas expansion thrust) [Newtons]
    ISP: float = 175.0                                                          # The estimated ISP of the engine [seconds]
    OF_ratio: float = 1.0                                                       # The target ratio of oxygen to fuel combustion in the engine [dimensionless] 
    total_mass_flow_rate: float = 3.82 * c.LBM2KG                               # The targeted mass flow rate through the engine [kilograms/second]
    oxidizer_mass_flow_rate: float = total_mass_flow_rate/(1 + (1/OF_ratio))    # The targeted mass flow rate for oxidizer through the engine [kilograms/second]
    core_fuel_mass_flow_rate: float = total_mass_flow_rate/(OF_ratio + 1)       # The targeted mass flow rate for fuel through the engine [kilograms/second]
    film_fuel_mass_flow_rate: float = 0.3 * core_fuel_mass_flow_rate            # The targeted film mass flow rate for fuel [kilograms/second]
    burn_time: float = 2.09                                                     # The estimated burn time of the engine [seconds]
    contraction_ratio: float = 7.0                                              # The target ratio of chamber area to throat area [dimensionless]
    exit_pressure: float = 15.0 * c.PSI2PA                                      # The target exit pressure of the exhaust gas [Pascals]
    # combustion_temperature: float = 2170                                      # The estimated combustion temperature [Kelvin]
    chamber_outer_diameter: float = 6.0 * c.IN2M                                # The design combustion chamber diameter [meters]
    chamber_inner_diameter: float = 4.9 * c.IN2M                                # The design combustion chamber diameter [meters]
    chamber_throat_diameter: float = 1.852 * c.IN2M                             # The design throat diameter [meters]
    
    # Tank Parameters
    # FYI the sizing of the tanks accounted for tank ullages and propellant residuals, so (burn_time * mass_flow_rate) will not equal total_propellant_mass.
    tank_pressure: float = 425 * c.PSI2PA     # The estimated required tank pressure to sustain the combustion pressure in the engine [Pascals]
    tank_outer_diameter: float = tube_outer_diameter  # Outer diameter of both tanks of the rocket
    tank_inner_diameter: float = tube_inner_diameter  # Inner diameter of both tanks of the rocket
    tank_wall_thickness: float = (tube_outer_diameter - tube_inner_diameter)/2  # Inner diameter of both tanks of the rocket
    
    fuel_tank_length: float = 6 * c.IN2M      # The length of the fuel tank that needs to be filled with fuel (the actual tank may be longer) [meters]
    fuel_tank_usable_volume: float = 2.55 * c.L2M3     # The required loaded volume of fuel needed for the burn time [meter^3]
    fuel_total_mass: float = 4.42 * c.LBM2KG     # The required loaded mass of fuel needed for the burn time [kilograms]
    fuel_used_mass: float = burn_time * core_fuel_mass_flow_rate 
    fuel_residual_mass: float = fuel_total_mass - fuel_used_mass

    oxidizer_tank_length: float = 4.56 * c.IN2M   # The length of the oxidizer tank that needs to be filled with oxidizer (the actual tank may be longer) [meters]
    oxidizer_tank_usable_volume: float = 1.94 * c.L2M3   # The required loaded volume of oxidizer needed for the burn time [meter^3]
    oxidizer_total_mass: float = 4.42 * c.LBM2KG   # The required loaded mass of oxidizer needed for the burn time [kilograms]
    oxidizer_used_mass: float = burn_time * oxidizer_mass_flow_rate 
    oxidizer_residual_mass: float = oxidizer_total_mass - oxidizer_used_mass

    total_propellant_mass: float = fuel_total_mass + oxidizer_total_mass # The total mass of propellant put on the vehicle [kilograms]
    total_used_propellant_mass: float = fuel_used_mass + oxidizer_used_mass # The total mass of propellant needed for the burn time [kilograms]
    total_residual_mass: float = fuel_residual_mass + oxidizer_residual_mass # The total mass of propellant put on the vehicle [kilograms]
    
    # all from CAD
    engine_length =         10.179 * c.IN2M
    injector_length =       0.475 * c.IN2M
    lower_length =          20 * c.IN2M
    bulkhead_length =       1.22 * c.IN2M
    mid_length =            5 * c.IN2M

    upper_length =          25 * c.IN2M
    recovery_bay_length =   24 * c.IN2M
    nosecone_length =       15 * c.IN2M

    # COPV Parameters
    COPV_volume: float = 4.70 * c.L2M3 # [m^3] 
    COPV_starting_pressure: float = 4300 * c.PSI2PA # [Pa]

    root_chord: float = 16 * c.IN2M
    tip_chord: float = 6 * c.IN2M
    sweep_length: float = root_chord - tip_chord
    wingspan: float = 7.6 * c.IN2M
    number_of_fins: float = 3 # [-]
    fin_top: float = None

    # fin_top = vehicle.lower_fuel_bulkhead.bottom_distance_from_aft
    # noseconeToFin = total_length - fin_top
        
    # 1-DoF Results:
    one_DoF_off_the_rail_TWR: float = 7.43                # The target thrust-to-weight ratio of the rocket off the launch rail [dimensionless]
    one_DoF_off_the_rail_acceleration: float = 6.43       # The target acceleration of the rocket off the launch rail [standard gravities]
    one_DoF_off_the_rail_velocity: float = 27.64          # The target velocity of the rocket off the launch rail [meters/second]

    one_DoF_max_acceleration: float = 6.96 # (upwards!)   # The maximum acceleration of the rocket during flight [standard gravities]
    one_DoF_max_mach: float = 0.395                       # The maximum speed of the rocket during flight [Mach (speed of sound of air)]
    one_DoF_max_velocity: float = one_DoF_max_mach * 343  # The maximum speed of the rocket during flight [meters/second]
    one_DoF_total_impulse: float = 6340                   # The total impulse of the rocket over the duration of flight [newton seconds]

    one_DoF_estimated_apogee: float = 2690 * c.FT2M       # The estimated 1-DoF altitude [meters]

    # 6-DoF results:
    # six_DoF_off_the_rail_TWR: float = ?                          # The target thrust-to-weight ratio of the rocket off the launch rail [dimensionless]
    six_DoF_off_the_rail_acceleration: float = 58.264 / c.GRAVITY  # The target acceleration of the rocket off the launch rail [standard gravities]
    six_DoF_off_the_rail_velocity: float = 28.70                   # The target velocity of the rocket off the launch rail [meters/second]

    # six_DoF_max_acceleration: float = ? # (upwards!)             # The maximum acceleration of the rocket during flight [standard gravities]
    six_DoF_max_mach: float = 0.368                                # The maximum speed of the rocket during flight [Mach (speed of sound of air)]
    six_DoF_max_velocity: float = six_DoF_max_mach * 343           # The maximum speed of the rocket during flight [meters/second]
    # six_DoF_total_impulse: float = ?                             # The total impulse of the rocket over the duration of flight [newton seconds]
    
    six_DoF_estimated_apogee: float = 747.60                       # The estimated 6-DoF altitude [meters]
    
    # later-calculated values that need to be here so it can still be added for a frozen data class
    wet_mass: float = None
    dry_mass: float = None
    total_length: float = None
    dry_COM_location_from_bottom: float = None
    dry_COM_location_from_top: float = None
    wet_COM_location_from_bottom: float = None
    wet_COM_location_from_top: float = None

    # internal freeze flag
    _frozen: bool = field(default=False, init=False, repr=False)

    def __post_init__(self):
        object.__setattr__(self, "_frozen", False)

    def freeze(self):
        object.__setattr__(self, "_frozen", True)

    def __setattr__(self, name, value):
        if getattr(self, "_frozen", False) and name != "_frozen":
            raise AttributeError("The vehicle parameters are frozen, you cannot change values, ask David how to change these values")
        super().__setattr__(name, value)

parameters = VehicleParameters()
  
# Center of mass !
# - 5:53 AM, 10/25/2025


# TO ADD:
# - changing mass over time
# - how do i add structures?
#     - strut mass


def CalcCylinderVolume(diameter, length):
    radius = diameter/2
    volume = np.pi * (radius**2) * length # off the dome!
    return volume

def CalcTubeVolume(OD, ID, length):
    volume = CalcCylinderVolume(OD, length) - CalcCylinderVolume(ID, length)
    return volume


propellant_tank_outer_diameter = parameters.tube_outer_diameter
propellant_tank_inner_diameter = parameters.tube_inner_diameter
panels_outer_diameter = parameters.tube_outer_diameter

panel_type = "foil"

if panel_type == "tube":
    panels_inner_diameter = parameters.tube_inner_diameter
elif panel_type == "foil":
    panels_thickness = 0.020 * c.IN2M
    panels_inner_diameter = panels_outer_diameter - (2 * panels_thickness)

# engine_wall_thickness = 0.25 * c.IN2M
# engine_OD = 6 * c.IN2M
# engine_ID = engine_OD - (2 * engine_wall_thickness)
# engine_mass = c.DENSITY_SS316 * CalcTubeVolume(engine_OD, engine_ID, engine_length)
engine_mass = 11.66 * c.LBM2KG # [lbm] measured CAD value

# injector_mass = c.DENSITY_SS316 * CalcCylinderVolume(propellant_tank_outer_diameter, injector_length)
injector_mass = 4.69 * c.LBM2KG # [lbm] measured CAD value

number_of_fins = 3
mass_per_fin = 1.614 * c.LBM2KG # [lbm]
total_fin_mass = number_of_fins * mass_per_fin # [kg]

number_of_fin_can_struts = 3
mass_per_fin_can_strut = 0.7257 * c.LBM2KG # [lbm]
total_lower_strut_mass = number_of_fin_can_struts * mass_per_fin_can_strut # [kg]

lower_panels_mass = c.DENSITY_AL * CalcTubeVolume(panels_outer_diameter, panels_inner_diameter, parameters.lower_length)

valves_mass = 2 * 3.26 * c.LBM2KG # fuel and ox 3/4 inch valve https://habonim.com/wp-content/uploads/2020/08/C47-BD_C47__2023_VO4_28-06-23.pdf

lower_mass = valves_mass + lower_panels_mass + total_fin_mass + total_lower_strut_mass

use_bulkhead_mass_estimate = False
if use_bulkhead_mass_estimate == True:
    bulkhead_wall_thickness = 0.25 * c.IN2M
    bulkhead_top_thickness = 0.76 * c.IN2M

    bulkhead_mass =  c.DENSITY_AL * (
        (CalcCylinderVolume(propellant_tank_outer_diameter, parameters.bulkhead_length) - 
        CalcCylinderVolume(propellant_tank_outer_diameter - (2 * bulkhead_wall_thickness), parameters.bulkhead_length - bulkhead_top_thickness))
    )
else:
    bulkhead_mass = 1.971 * c.LBM2KG # [lbm] measured CAD value

fuel_tank_wall_mass = c.DENSITY_AL * CalcTubeVolume(propellant_tank_outer_diameter, propellant_tank_inner_diameter, parameters.fuel_tank_length)
fuel_tank_dry_mass = fuel_tank_wall_mass + parameters.fuel_residual_mass
fuel_tank_wet_mass = fuel_tank_dry_mass + parameters.fuel_used_mass

# total_mid_strut_mass
mid_panels_mass = c.DENSITY_AL * CalcTubeVolume(panels_outer_diameter, panels_inner_diameter, parameters.mid_length)
mid_mass = mid_panels_mass

oxidizer_tank_wall_mass = c.DENSITY_AL * CalcTubeVolume(propellant_tank_outer_diameter, propellant_tank_inner_diameter, parameters.oxidizer_tank_length)
oxidizer_tank_dry_mass = oxidizer_tank_wall_mass + parameters.oxidizer_residual_mass
oxidizer_tank_wet_mass = oxidizer_tank_dry_mass + parameters.oxidizer_used_mass

regulator_mass = 1.200 # regulator https://valvesandregulators.aquaenvironment.com/item/high-flow-reducing-regulators-2/873-d-high-flow-dome-loaded-reducing-regulators/item-1659
copv_mass = 2.9 # [kg]
upper_airframe_tube_mass = c.DENSITY_AL * CalcTubeVolume(parameters.tube_outer_diameter, parameters.tube_inner_diameter, parameters.upper_length)
upper_mass = regulator_mass + copv_mass + upper_airframe_tube_mass

recovery_bay_airframe_tube_mass = c.DENSITY_AL * CalcTubeVolume(panels_outer_diameter, panels_inner_diameter, parameters.recovery_bay_length)
parachute_mass = 2.25 * c.LBM2KG  # [kg] https://shop.fruitychutes.com/collections/parachutes/products/iris-ultra-144-compact-chute-114lbs-20fps-64lbs-15fps
recovery_bay_mass = recovery_bay_airframe_tube_mass + parachute_mass

tungsten_cube = 10 * c.LBM2KG

nosecone_mass = c.DENSITY_AL * CalcTubeVolume(panels_outer_diameter, panels_inner_diameter, parameters.nosecone_length) + tungsten_cube

# structures = 15 * c.LBM2KG # structures ! funny

@dataclass(frozen=True)
class MassComponent:
    name : str                # name of the component [string]
    mass: float               # [kilograms]
    bottom_distance_from_aft: float   # distance from the bottom of the rocket to the bottom of the mass component [meters]
    length: float              # [meters]
    top_distance_from_aft: float = None

    def __post_init__(self):
        object.__setattr__(
            self,
            "top_distance_from_aft",
            self.bottom_distance_from_aft + self.length
        )

@dataclass(frozen=True)
class MassDistribution:
    components: list[MassComponent]

    def __iter__(self):
        return iter(self.components)

engine =                  MassComponent(name = "engine",                      mass = engine_mass,            bottom_distance_from_aft = 0,                                                   length = parameters.engine_length)
injector =                MassComponent(name = "injector",                    mass = injector_mass,          bottom_distance_from_aft = engine.top_distance_from_aft,                        length = parameters.injector_length)
lower =                   MassComponent(name = "lower",                       mass = lower_mass,             bottom_distance_from_aft = injector.top_distance_from_aft,                      length = parameters.lower_length)

lower_fuel_bulkhead =     MassComponent(name = "lower_fuel_bulkhead",         mass = bulkhead_mass,          bottom_distance_from_aft = lower.top_distance_from_aft,                         length = parameters.bulkhead_length)
wet_fuel_tank =           MassComponent(name = "wet_fuel_tank",               mass = fuel_tank_wet_mass,     bottom_distance_from_aft = lower_fuel_bulkhead.top_distance_from_aft,           length = parameters.fuel_tank_length)
dry_fuel_tank =           MassComponent(name = "dry_fuel_tank",               mass = fuel_tank_dry_mass,     bottom_distance_from_aft = lower_fuel_bulkhead.top_distance_from_aft,           length = parameters.fuel_tank_length)
upper_fuel_bulkhead =     MassComponent(name = "upper_fuel_bulkhead",         mass = bulkhead_mass,          bottom_distance_from_aft = wet_fuel_tank.top_distance_from_aft,                 length = parameters.bulkhead_length)

mid =                     MassComponent(name = "mid",                         mass = mid_mass,               bottom_distance_from_aft = upper_fuel_bulkhead.top_distance_from_aft,           length = parameters.mid_length)

lower_oxidizer_bulkhead = MassComponent(name = "lower_oxidizer_bulkhead",     mass = bulkhead_mass,          bottom_distance_from_aft = mid.top_distance_from_aft,                           length = parameters.bulkhead_length)
wet_oxidizer_tank =       MassComponent(name = "wet_oxidizer_tank",           mass = oxidizer_tank_wet_mass, bottom_distance_from_aft = lower_oxidizer_bulkhead.top_distance_from_aft,       length = parameters.oxidizer_tank_length)
dry_oxidizer_tank =       MassComponent(name = "dry_oxidizer_tank",           mass = oxidizer_tank_dry_mass, bottom_distance_from_aft = lower_oxidizer_bulkhead.top_distance_from_aft,       length = parameters.oxidizer_tank_length)
upper_oxidizer_bulkhead = MassComponent(name = "upper_oxidizer_bulkhead",     mass = bulkhead_mass,          bottom_distance_from_aft = wet_oxidizer_tank.top_distance_from_aft,             length = parameters.bulkhead_length)

upper =                   MassComponent(name = "upper",                       mass = upper_mass,             bottom_distance_from_aft = upper_oxidizer_bulkhead.top_distance_from_aft,       length = parameters.upper_length)
recovery_bay =            MassComponent(name = "recovery_bay",                mass = recovery_bay_mass,      bottom_distance_from_aft = upper.top_distance_from_aft,                         length = parameters.recovery_bay_length)

nosecone =                MassComponent(name = "nosecone",                    mass = nosecone_mass,          bottom_distance_from_aft = recovery_bay.top_distance_from_aft,                  length = parameters.nosecone_length)


wet_mass_distribution = MassDistribution(components=
    [
    engine,
    injector,
    lower,
    
    lower_fuel_bulkhead,
    wet_fuel_tank,
    upper_fuel_bulkhead,
    
    mid,
    
    lower_oxidizer_bulkhead,
    wet_oxidizer_tank,
    upper_oxidizer_bulkhead,
    
    upper,
    
    recovery_bay,
    nosecone,
    ]
)

dry_mass_distribution = MassDistribution(components=
    [
    engine,
    injector,
    lower,
    
    lower_fuel_bulkhead,
    dry_fuel_tank,
    upper_fuel_bulkhead,
    
    mid,
    
    lower_oxidizer_bulkhead,
    dry_oxidizer_tank,
    upper_oxidizer_bulkhead,
    
    upper,
    
    recovery_bay,
    nosecone,
    ]
)



def calcCG(linear_density_array, length_along_rocket_linspace):
    '''
    linear_density_array: Array of linear density as a function of length [kg / m]
    length_along_rocket_linspace: Array of length along rocket [m]
    cg: Location of center of gravity of rocket from aft [m]
    '''
    dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]
    totalMass = np.sum(linear_density_array * dx)
    # print(totalMass / LB2KG)
    lengths = np.array(length_along_rocket_linspace)
    masses = np.array(linear_density_array * dx)
    moments = np.sum(lengths * masses)
    cg = moments / totalMass
    return cg

rocket_dict_wet = {} # Rocket dictionary for wet mass
for item in wet_mass_distribution.components:
    rocket_dict_wet[item.name] = {
        "mass": item.mass,
        "bottom_distance_from_aft": item.bottom_distance_from_aft,
        "length": item.length
    }

rocket_dict_dry = {} # Rocket dictionary for dry mass
for item in dry_mass_distribution.components:
    rocket_dict_dry[item.name] = {
        "mass": item.mass,
        "bottom_distance_from_aft": item.bottom_distance_from_aft,
        "length": item.length
    }


num_points = 2000
length_along_rocket_linspace = np.linspace(wet_mass_distribution.components[0].bottom_distance_from_aft, wet_mass_distribution.components[-1].top_distance_from_aft, num_points)

wet_linear_density_array = np.zeros(num_points)
dry_linear_density_array = np.zeros(num_points)

for wet_component, dry_component in zip(wet_mass_distribution.components, dry_mass_distribution.components, strict=True):
    
    wet_linear_density = (wet_component.mass / wet_component.length) # The average mass in the length of the component
    dry_linear_density = (dry_component.mass / dry_component.length) # The average mass in the length of the component
    
    # should be the same for wet and dry so i just used wet
    for index, length_along_rocket in enumerate(length_along_rocket_linspace):
        above_component_bottom = length_along_rocket >= wet_component.bottom_distance_from_aft
        below_component_top = length_along_rocket <= (wet_component.bottom_distance_from_aft + wet_component.length)
        
        if (above_component_bottom and below_component_top):
            wet_linear_density_array[index] += wet_linear_density
            dry_linear_density_array[index] += dry_linear_density

vehicle_wet_mass = sum(component.mass for component in wet_mass_distribution)
vehicle_dry_mass = sum(component.mass for component in dry_mass_distribution)

# check that last length value is the longest value (would be true if code works as intended)
if max(length_along_rocket_linspace) != length_along_rocket_linspace[-1]:
    raise ValueError("somethings wrong...")

rocket_length = max(length_along_rocket_linspace)

parameters.fin_top = injector.top_distance_from_aft + parameters.root_chord # assume fins start right after injector

parameters.wet_mass = vehicle_wet_mass # The estimated dry mass of the rocket [kilograms]
parameters.dry_mass = vehicle_dry_mass # The estimated dry mass of the rocket [kilograms]
parameters.total_length = rocket_length # The estimated length of the rocket [meters]

parameters.wet_COM_location_from_bottom = calcCG(wet_linear_density_array, length_along_rocket_linspace)
parameters.wet_COM_location_from_top = rocket_length - parameters.wet_COM_location_from_bottom

parameters.dry_COM_location_from_bottom = calcCG(dry_linear_density_array, length_along_rocket_linspace)
parameters.dry_COM_location_from_top = rocket_length - parameters.dry_COM_location_from_bottom

parameters.freeze()
# parameters.wet_mass = 9999999999999999999999999999999999999999 # if you wanna test that freezing works




########## paused until i figure ts out ##########
# repository_root_path = vehicle_parameters_functions.Get_Repository_Root_Path()

# PSPL_ROCKET_A_file_path = repository_root_path / f"vehicle_parameters.csv"

# export_path_list = [PSPL_ROCKET_A_file_path]

# try:
#     Six_DoF_csv_file_path = (
#         repository_root_path
#         / ".."          # one directory up from PSPL_Rocket_A
#         / "PSPL-6DOF"
#         / "TheSixDoF"
#         / "Inputs"
#         / "Saved Rockets"
#         / "FUCK_MATLAB"
#         / "vehicle_parameters.csv"
#     ).resolve()

#     # write to CSV for 6DOF to read since 6DOF is in matlabese
#     Six_DoF_csv_file_path.parent.mkdir(parents=True, exist_ok=True)
#     export_path_list.append(Six_DoF_csv_file_path)
# except:
#     pass




# for export_file_path in export_path_list:
#     vehicle_parameters_functions.ConvertObjectToCSV(wet_mass_distribution, "vehicle_wet_mass_distribution")

#     is_equal = vehicle_parameters_functions.Determine_if_CSV_Files_are_Equal(
#         "wet_mass_distribution.csv",
#         "wet_mass_distribution_new.csv",
#     )

#     if not is_equal:
#         pass







if __name__ == "__main__":
    
    print(f"Vehicle Wet Mass: {vehicle_wet_mass * c.KG2LBM:.2f} lbm, {vehicle_wet_mass:.2f} kg")
    print(f"Vehicle Dry Mass: {vehicle_dry_mass * c.KG2LBM:.2f} lbm, {vehicle_dry_mass:.2f} kg")
    
    panels_mass = lower_panels_mass + mid_panels_mass
    # print(f"panels mass: {panels_mass * c.KG2LBM:.2f} lbm")

    
    print(f"\nRocket Length: {rocket_length * c.M2IN:.2f} in, {rocket_length * c.M2FT:.2f} ft")
    print(f"Rocket Length: {rocket_length:.2f} m\n")
    
    print(f"Wet CoM location distance from bottom: {parameters.wet_COM_location_from_bottom * c.M2IN:.2f} in, {parameters.wet_COM_location_from_bottom:.3f} m")
    print(f"Dry CoM location distance from bottom: {parameters.dry_COM_location_from_bottom * c.M2IN:.2f} in, {parameters.dry_COM_location_from_bottom:.3f} m")
    
    print(f"\nWet CoM location distance from top:    {parameters.wet_COM_location_from_top * c.M2IN:.2f} in, {parameters.wet_COM_location_from_top:.3f} m")
    print(f"Dry CoM location distance from top:    {parameters.dry_COM_location_from_top * c.M2IN:.2f} in, {parameters.dry_COM_location_from_top:.3f} m")
    
    plt.plot(length_along_rocket_linspace * c.M2FT, (wet_linear_density_array * (c.KG2LBM / c.M2FT)))
    plt.vlines(parameters.wet_COM_location_from_bottom * c.M2FT, min(wet_linear_density_array * (c.KG2LBM / c.M2FT)), max(wet_linear_density_array * (c.KG2LBM / c.M2FT)), color="red", linestyles="dotted", label="Center of Gravity")
    plt.legend()
    
    plt.xlabel("Length from Bottom [feet]")
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    
    plt.ylabel("Wet Mass Density [lbs/feet^3]")
    plt.title("Rocket Wet Mass Distribution")


    print_components = False
    
    if print_components == True:
        for component in wet_mass_distribution.components:
            print(f"{component.name}:")
            
            # imperial
            print(f"\tlength: {component.length * c.M2IN:.2f} in")
            print(f"\tmass: {component.mass * c.KG2LBM:.2f} lbm")
            print(f"\tdistance from top: {(rocket_length - (component.bottom_distance_from_aft + (component.length/2))) * c.M2IN:.2f} in")
            
            # # metric
            
            # component.StartAfter()
            # airframe_length = rocket_length - nosecone.length
            # six_dof_bottom = airframe_length - component.StartAfter()
            
            # six_dof_middle = six_dof_bottom + (component.length/2)
            # print(f"\t6dof x bottom: {six_dof_bottom + component.length:.2f} m")
            # print(f"\t6dof x middle: {six_dof_middle:.2f} m")
            
            # print(f"\tlength: {component.length:.2f} m")
            # print(f"\tmass: {component.mass:.2f} kg")
            # print(f"\tdistance from top: {(rocket_length - (component.bottom_distance_from_aft + (component.length/2))):.2f} m")
            
    # plt.show()