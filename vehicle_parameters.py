# This code reflects the vehicle parameters as stated on the vehicle parameters page: https://purdue-space-program.atlassian.net/wiki/x/Aw5qYw
# Import this file into your script as a way to easily and reliably pull parameters. 
# Author: David Gustafsson

from dataclasses import dataclass, fields
from typing import Iterator
import matplotlib.pyplot as plt
import numpy as np

import constants as c



@dataclass(frozen=True)
class VehicleParameters:

    fuel_name: str = "isopropyl alcohol" 
    oxidizer_name: str = "liquid oxygen"
    tube_outer_diameter: float = 6.0 * c.IN2M         # Outer diameter of tube used in some sections of the rocket
    tube_inner_diameter: float = 5.75 * c.IN2M        # Inner diameter of tube used in some sections of the rocket
    
    
    chamber_pressure: float = 250 * c.PSI2PA                                    # The target combustion pressure in the engine [Pascals]
    jet_thrust: float = 668.0 * c.LBF2N                                         # The targeted engine thrust (not accounting for exhaust gas expansion thrust) [Newtons]
    ISP: float = 175.0                                                          # The estimated ISP of the engine [seconds]
    OF_ratio: float = 1.0                                                       # The target ratio of oxygen to fuel combustion in the engine [dimensionless] 
    total_mass_flow_rate: float = 3.82 * c.LBM2KG                               # The targeted mass flow rate through the engine [kilograms/second]
    oxidizer_mass_flow_rate: float = total_mass_flow_rate/(1 + (1/OF_ratio))     # The targeted mass flow rate for oxidizer through the engine [kilograms/second]
    fuel_mass_flow_rate: float = total_mass_flow_rate/(OF_ratio + 1)            # The targeted mass flow rate for fuel through the engine [kilograms/second]
    burn_time: float = 2.09                                                     # The estimated burn time of the engine [seconds]
    contraction_ratio: float = 7.0                                              # The target ratio of chamber area to throat area [dimensionless]
    exit_pressure: float = 15.0 * c.PSI2PA                                      # The target exit pressure of the exhaust gas [Pascals]
    # combustion_temperature: float = 2170                                      # The estimated combustion temperature [Kelvin]
    chamber_outer_diameter: float = 6.0 * c.IN2M                                # The design combustion chamber diameter [meters]
    chamber_inner_diameter: float = 4.9 * c.IN2M                                # The design combustion chamber diameter [meters]
    chamber_throat_diameter: float = 1.852 * c.IN2M                             # The design throat diameter [meters]
    
    # FYI the sizing of the tanks accounted for tank ullages and propellant residuals, so (burn_time * mass_flow_rate) will not equal total_propellant_mass.

    tank_outer_diameter: float = tube_outer_diameter  # Outer diameter of both tanks of the rocket
    tank_inner_diameter: float = tube_inner_diameter  # Inner diameter of both tanks of the rocket
    tank_wall_thickness: float = (tube_outer_diameter - tube_inner_diameter)/2  # Inner diameter of both tanks of the rocket
    tank_pressure: float = 416.67 * c.PSI2PA     # The estimated required tank pressure to sustain the combustion pressure in the engine [Pascals]
    
    fuel_tank_length: float = 6 * c.IN2M      # The length of the fuel tank that needs to be filled with fuel (the actual tank may be longer) [meters]
    fuel_tank_usable_volume: float = 2.55 * c.L2M3     # The required loaded volume of fuel needed for the burn time [meter^3]
    fuel_total_mass: float = 4.42 * c.LBM2KG     # The required loaded mass of fuel needed for the burn time [kilograms]

    oxidizer_tank_length: float = 4.56 * c.IN2M   # The length of the oxidizer tank that needs to be filled with oxidizer (the actual tank may be longer) [meters]
    oxidizer_tank_usable_volume: float = 1.94 * c.L2M3   # The required loaded volume of oxidizer needed for the burn time [meter^3]
    oxidizer_total_mass: float = 4.42 * c.LBM2KG   # The required loaded mass of oxidizer needed for the burn time [kilograms]

    total_propellant_mass: float = fuel_total_mass + oxidizer_total_mass # (4.42 + 4.42) * c.LB2KG # The total mass of propellant needed for the burn time [kilograms]

    total_length: float = 7.5 * c.FT2M            # The estimated length of the rocket [meter]
    wet_mass: float = 93.3 * c.LBM2KG             # The estimated dry mass of the rocket [kilograms]
    dry_mass: float = 84.5 * c.LBM2KG             # The estimated dry mass of the rocket [kilograms]
    estimated_apogee: float = 2690 * c.FT2M       # The estimated 1-DOF altitude [meters]
    off_the_rail_TWR: float = 7.43                # The target thrust-to-weight ratio of the rocket off the launch rail [dimensionless]
    off_the_rail_acceleration: float = 6.43       # The target acceleration of the rocket off the launch rail [standard gravity]
    off_the_rail_velocity: float = 27.64          # The target velocity of the rocket off the launch rail [meters/second]
    max_acceleration: float = 6.96 # (upwards!)   # The maximum acceleration of the rocket during flight [standard gravity]
    max_mach: float = 0.395                       # The maximum speed of the rocket during flight [Mach (speed of sound of air)]
    max_velocity: float = max_mach * 343          # The maximum speed of the rocket during flight [meters/second]
    total_impulse: float = 6340                   # The total impulse of the rocket over the duration of flight [newton seconds]
    
parameters = VehicleParameters()


# vehicle parameters that did not come from the Rocket A sizing code
@dataclass(frozen=True)
class CalculatedVehicleParameters:
    COPV_volume: float = 4.70 * c.L2M3 # [m^3] 
    COPV_starting_pressure: float =  4300 * c.PSI2PA # [Pa]
    

calculated_parameters = CalculatedVehicleParameters()

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
    
engine_length =         10.179 * c.IN2M
injector_length =       0.475 * c.IN2M
lower_length =          12 * c.IN2M
bulkhead_length =       2 * c.IN2M
fuel_tank_length =      parameters.fuel_tank_length
mid_length =            5 * c.IN2M
oxidizer_tank_length =  parameters.oxidizer_tank_length

upper_length =          12 * c.IN2M
helium_bay_length =     16 * c.IN2M

avionics_bay_length =   3 * c.IN2M
recovery_bay_length =   6 * c.IN2M
nosecone_length =       12 * c.IN2M

propellant_tank_outer_diameter = parameters.tube_outer_diameter
propellant_tank_inner_diameter = parameters.tube_inner_diameter
panels_outer_diameter = parameters.tube_outer_diameter
panels_inner_diameter = parameters.tube_inner_diameter

engine_wall_thickness = 0.25 * c.IN2M
# engine_ID = propellant_tank_outer_diameter - (2 * engine_wall_thickness)
engine_OD = 6 * c.IN2M
engine_ID = engine_OD - (2 * engine_wall_thickness)
engine_mass = c.DENSITY_SS316 * CalcTubeVolume(engine_OD, engine_ID, engine_length)

injector_mass = c.DENSITY_SS316 * CalcCylinderVolume(propellant_tank_outer_diameter, injector_length)

number_of_fins = 3
fin_mass = number_of_fins * 3.51 * c.LBM2KG # [kg]
valves_mass = 2 * 3.26 * c.LBM2KG # fuel and ox 3/4 inch valve https://habonim.com/wp-content/uploads/2020/08/C47-BD_C47__2023_VO4_28-06-23.pdf
lower_panels_mass = c.DENSITY_AL * CalcTubeVolume(panels_outer_diameter, panels_inner_diameter, lower_length)
lower_mass = valves_mass + lower_panels_mass + fin_mass

bulkhead_wall_thickness = 0.25 * c.IN2M
bulkhead_top_thickness = 0.76 * c.IN2M

bulkhead_mass =  c.DENSITY_AL * (
    (CalcCylinderVolume(propellant_tank_outer_diameter, bulkhead_length) - 
    CalcCylinderVolume(propellant_tank_outer_diameter - (2 * bulkhead_wall_thickness), bulkhead_length - bulkhead_top_thickness))
)

fuel_tank_wall_mass = c.DENSITY_AL * CalcTubeVolume(propellant_tank_outer_diameter, propellant_tank_inner_diameter, fuel_tank_length)
fuel_tank_wet_mass = fuel_tank_wall_mass + parameters.fuel_total_mass


mid_panels_mass = c.DENSITY_AL * CalcTubeVolume(panels_outer_diameter, panels_inner_diameter, mid_length)
mid_mass = mid_panels_mass

oxidizer_tank_wall_mass = c.DENSITY_AL * CalcTubeVolume(propellant_tank_outer_diameter, propellant_tank_inner_diameter, oxidizer_tank_length)
oxidizer_tank_wet_mass = oxidizer_tank_wall_mass + parameters.oxidizer_total_mass

regulator_mass = 1.200 # regulator https://valvesandregulators.aquaenvironment.com/item/high-flow-reducing-regulators-2/873-d-high-flow-dome-loaded-reducing-regulators/item-1659
upper_panels_mass = c.DENSITY_AL * CalcTubeVolume(panels_outer_diameter, panels_inner_diameter, upper_length)
upper_mass = regulator_mass + upper_panels_mass

copv_mass = 2.9 
helium_bay_panels_mass = c.DENSITY_AL * CalcTubeVolume(panels_outer_diameter, panels_inner_diameter, helium_bay_length)
helium_bay_mass = copv_mass + helium_bay_panels_mass

avionics_bay_panels_mass = c.DENSITY_AL * CalcTubeVolume(panels_outer_diameter, panels_inner_diameter, avionics_bay_length)
avionics_bay_mass = avionics_bay_panels_mass + (1 * c.LBM2KG) # avionics doesn't weigh anything...

recovery_bay_panels_mass = c.DENSITY_AL * CalcTubeVolume(panels_outer_diameter, panels_inner_diameter, recovery_bay_length)
parachute_mass = 8 * c.LBM2KG  # [kg] 1/3 cuz 1/3 of dry mass compared to --> https://github.com/Purdue-Space-Program/PSPL_Rocket_4_Sizing/blob/2b15e1dc508a56731056ff594a3c6b5afb639b4c/scripts/structures.py#L75
recovery_bay_mass = recovery_bay_panels_mass + parachute_mass

nose_cone_mass = c.DENSITY_AL * CalcTubeVolume(panels_outer_diameter, panels_inner_diameter, nosecone_length) # guess



structures = 15 * c.LBM2KG # structures !



@dataclass(frozen=True)
class MassComponent:
    name : str                # name of the component [string]
    mass: float               # [kilograms]
    bottom_distance_from_aft: float   # distance from the bottom of the rocket to the bottom of the mass component [meters]
    length: float              # [meters]
    
    def StartAfter(self): # A mass object that directly after another
        return(self.bottom_distance_from_aft + self.length)

@dataclass(frozen=True)
class MassDistribution:
    components: list[MassComponent]

    def __iter__(self):
        return iter(self.components)

    # def __iter__(self) -> Iterator[MassComponent]:
    #     """Automatically iterate over all MassComponent attributes."""
    #     for f in fields(self):
    #         yield getattr(self, f.name)


engine =                  MassComponent(name = 'engine',                      mass = engine_mass,            bottom_distance_from_aft = 0,                                                length = engine_length)
injector =                MassComponent(name = 'injector',                    mass = injector_mass,          bottom_distance_from_aft = engine.StartAfter(),                              length = injector_length)
lower =                   MassComponent(name = 'lower',                       mass = lower_mass,             bottom_distance_from_aft = injector.StartAfter(),                            length = lower_length)

lower_fuel_bulkhead =     MassComponent(name = 'lower_fuel_bulkhead',         mass = bulkhead_mass,          bottom_distance_from_aft = lower.StartAfter(),                               length = bulkhead_length)
fuel_tank =               MassComponent(name = 'fuel_tank',                   mass = fuel_tank_wet_mass,     bottom_distance_from_aft = lower_fuel_bulkhead.StartAfter(),                 length = parameters.fuel_tank_length)
upper_fuel_bulkhead =     MassComponent(name = 'upper_fuel_bulkhead',         mass = bulkhead_mass,          bottom_distance_from_aft = fuel_tank.StartAfter() - (bulkhead_length),       length = bulkhead_length)

mid =                     MassComponent(name = 'mid',                         mass = mid_mass,               bottom_distance_from_aft = upper_fuel_bulkhead.StartAfter(),                 length = mid_length)

lower_oxidizer_bulkhead = MassComponent(name = 'lower_oxidizer_bulkhead',     mass = bulkhead_mass,          bottom_distance_from_aft = mid.StartAfter(),                                 length = bulkhead_length)
oxidizer_tank =           MassComponent(name = 'oxidizer_tank',               mass = oxidizer_tank_wet_mass, bottom_distance_from_aft = lower_oxidizer_bulkhead.StartAfter(),             length = parameters.oxidizer_tank_length)
upper_oxidizer_bulkhead = MassComponent(name = 'upper_oxidizer_bulkhead',     mass = bulkhead_mass,          bottom_distance_from_aft = oxidizer_tank.StartAfter() - (bulkhead_length),   length = bulkhead_length)

upper =                   MassComponent(name = 'upper',                       mass = upper_mass,             bottom_distance_from_aft = upper_oxidizer_bulkhead.StartAfter(),             length = upper_length)
helium_bay =              MassComponent(name = 'helium_bay',                  mass = helium_bay_mass,        bottom_distance_from_aft = upper.StartAfter(),                               length = helium_bay_length)
avionics_bay =            MassComponent(name = 'avionics_bay',                mass = avionics_bay_mass,      bottom_distance_from_aft = helium_bay.StartAfter(),                          length = avionics_bay_length)
recovery_bay =            MassComponent(name = 'recovery_bay',                mass = recovery_bay_mass,      bottom_distance_from_aft = avionics_bay.StartAfter(),                        length = recovery_bay_length)

nosecone =                MassComponent(name = 'nosecone',                    mass = nose_cone_mass,         bottom_distance_from_aft = recovery_bay.StartAfter(),                        length=nosecone_length)


mass_distribution = MassDistribution(components=
    [
    engine,
    injector,
    lower,
    
    lower_fuel_bulkhead,
    fuel_tank,
    upper_fuel_bulkhead,
    
    lower_oxidizer_bulkhead,
    oxidizer_tank,
    upper_oxidizer_bulkhead,
    
    upper,
    helium_bay,
    
    avionics_bay,
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
for item in mass_distribution.components:
    rocket_dict_wet[item.name] = {
        "mass": item.mass,
        "bottom_distance_from_aft": item.bottom_distance_from_aft,
        "length": item.length
    }

rocket_dict_dry = {} # Rocket dictionary for dry mass
for item in rocket_dict_wet:
    if item == "fuel_tank":
        rocket_dict_dry[item] = {
            "mass": rocket_dict_wet[item]["mass"] - parameters.fuel_total_mass,
            "bottom_distance_from_aft": rocket_dict_wet[item]["bottom_distance_from_aft"],
            "length": rocket_dict_wet[item]["length"]
        }
    elif item == "oxidizer_tank":
        rocket_dict_dry[item] = {
            "mass": rocket_dict_wet[item]["mass"] - parameters.oxidizer_total_mass,
            "bottom_distance_from_aft": rocket_dict_wet[item]["bottom_distance_from_aft"],
            "length": rocket_dict_wet[item]["length"]
        }
    else:
        rocket_dict_dry[item] = rocket_dict_wet[item]
item_sum = 0
for item in rocket_dict_wet: item_sum += rocket_dict_wet[item]['mass']; print(f"{item}: {rocket_dict_wet[item]['mass']*c.KG2LBM} lbm")
print(f"Mass of rocket dict wet: {item_sum * c.KG2LBM} lbm")
for item in rocket_dict_wet: print(f"{item}: {rocket_dict_dry[item]['length']*c.M2FT} ft")

'''
print(rocket_dict_wet)
print(f"Mass of rocket dict wet: {sum(component['mass'] for component in rocket_dict_wet.values())} kg")
print(" ")
print(rocket_dict_dry)
print(f"Mass of rocket dict dry: {sum(component['mass'] for component in rocket_dict_dry.values())} kg")
print(" ")
print(rocket_dict_recovery)
print(f"Mass of rocket dict recovery: {sum(component['mass'] for component in rocket_dict_recovery.values())} kg")
'''
num_points = 500
length_along_rocket_linspace = np.linspace(mass_distribution.components[0].bottom_distance_from_aft, mass_distribution.components[-1].StartAfter(), num_points)

# x = np.linspace(0, nosecone.StartAfter(), num_points_per_component * np.size(mass_distribution))
# y = np.zeros(num_points_per_component * len(mass_distribution))

linear_density_array = np.zeros(num_points)

for component in mass_distribution.components:
    
    linear_density = (component.mass / component.length) # The average mass in the length of the component
    
    for index, length_along_rocket in enumerate(length_along_rocket_linspace):
        above_component_bottom = length_along_rocket >= component.bottom_distance_from_aft
        below_component_top = length_along_rocket <= (component.bottom_distance_from_aft + component.length)
        
        if (above_component_bottom and below_component_top):
            linear_density_array[index] += linear_density


if __name__ == "__main__":
    print(f"wet mass: {sum(component.mass for component in mass_distribution) * c.KG2LBM} lbm")
    print(f"engine mass: {engine.mass * c.KG2LBM:.2f} lbm")
    print(f"injector mass: {injector.mass * c.KG2LBM:.2f} lbm")

    panels_mass = lower_panels_mass + upper_panels_mass + helium_bay_panels_mass + avionics_bay_panels_mass + recovery_bay_panels_mass
    print(f"panels mass: {panels_mass * c.KG2LBM:.2f} lbm")
    print(f"fuel tank mass: {fuel_tank_wall_mass * c.KG2LBM:.2f} lbm")
    print(f"oxidizer tank mass: {oxidizer_tank_wall_mass * c.KG2LBM:.2f} lbm")
    print(f"nose cone mass: {nose_cone_mass * c.KG2LBM:.2f} lbm")
    print(f"upper mass: {upper_mass * c.KG2LBM:.2f} lbm")
    print(f"mid mass: {mid_mass * c.KG2LBM:.2f} lbm")
    print(f"lower mass: {lower_mass * c.KG2LBM:.2f} lbm")
    print(f"helium bay mass: {helium_bay_mass * c.KG2LBM:.2f} lbm")
    print(f"recovery bay mass: {recovery_bay_mass * c.KG2LBM:.2f} lbm")
    print(f"fuel mass: {parameters.fuel_total_mass * c.KG2LBM:.2f} lbm")


    plt.plot(length_along_rocket_linspace * c.M2FT, (linear_density_array * (c.KG2LBM / c.M2FT))    )
    
    rocket_length = max(length_along_rocket_linspace)
    print(f"\nRocket Length: {rocket_length * c.M2IN:.2f} in, {rocket_length * c.M2FT:.2f} ft")
    
    COG_location = calcCG(linear_density_array, length_along_rocket_linspace)
    print(f"COM location distance from bottom: {COG_location * c.M2IN:.2f} in")
    print(f"COM location distance from top: {(rocket_length - COG_location) * c.M2IN:.2f} in")
    
    plt.vlines(COG_location * c.M2FT, min(linear_density_array * (c.KG2LBM / c.M2FT)), max(linear_density_array * (c.KG2LBM / c.M2FT)), color="red", linestyles="dotted", label="Center of Gravity")
    plt.legend()
    
    plt.xlabel("length [feet]")
    plt.ylabel("mass density [lbs/feet]")
    plt.show()

    for component in mass_distribution.components:
        print(f"{component.name}: {component.mass:.2f} kg {component.length:.2f} m long")
        print(f"\t{component.mass:.2f} kg")
        print(f"\t{(rocket_length - (component.bottom_distance_from_aft + (component.length/2))):.2f} m from nose")
    print(f"Fuel mass: {parameters.fuel_total_mass:.2f} kg")
    print(f"fuel tank mass: {fuel_tank_wall_mass+lower_fuel_bulkhead.mass+upper_fuel_bulkhead.mass:.2f} kg")
    print(f"total fuel tank length: {fuel_tank.length + lower_fuel_bulkhead.length + upper_fuel_bulkhead.length:.2f} m")
    print(f"Oxidizer mass: {parameters.oxidizer_total_mass:.2f} kg")
    print(f"oxidizer tank mass: {oxidizer_tank_wall_mass+lower_oxidizer_bulkhead.mass+upper_oxidizer_bulkhead.mass:.2f} kg")
    print(f"total oxidizer tank length: {oxidizer_tank.length + lower_oxidizer_bulkhead.length + upper_oxidizer_bulkhead.length:.2f} m")

    
