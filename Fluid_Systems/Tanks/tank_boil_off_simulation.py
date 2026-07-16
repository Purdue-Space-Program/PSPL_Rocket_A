import sys
import os
import glob
import importlib

import sys,pathlib,collections,importlib.abc; r=next(p for p in pathlib.Path(__file__).resolve().parents if p.name=="PSPL_Rocket_A");   m=collections.defaultdict(list); [m[f.stem.casefold()].append(f) for f in r.rglob("*.py") if f.stem.isidentifier() and f.name != "__init__.py" and not (set(f.relative_to(r).parts)&{".git",".venv","__pycache__","build","dist"})]; dup={k:v for k,v in m.items() if   len(v)>1}; sys.meta_path.insert(0,type("AmbiguousBareImportBlocker",(importlib.abc.MetaPathFinder,),{"find_spec":lambda   self,fullname,path=None,target=None: (_ for _ in ()).throw(ImportError(f"The import {fullname!r} could refer to any following packages: "+" ".join("\n\t" + str(p.relative_to(r)) for p in dup[fullname.casefold()])+f"\n\nSpecify which package it is by using the folder.\nFor example:\n\t'import SFD.{fullname}'")) if "." not in fullname and fullname.casefold() in dup else None})()); sys.path.insert(0,str(r)) if str(r) not   in sys.path else None; [sys.path.append(str(v[0].parent)) for k,v in m.items() if k not in dup and str(v[0].parent) not in sys.path]

import vehicle_parameters_functions
import vehicle_parameters
import vehicle_main
import print_filter

import constants as c
import numpy as np
import CoolProp
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt
from dataclasses import dataclass, field


def Simulate_Boil_Off(parameters):
    
    def Calculate_Circle_Area(diameter):
        circle_radius = diameter/2
        circle_area = np.pi*(circle_radius**2)
        return (circle_area)
    
    def Calculate_Air_Resistance(outer_surface_area):
        air_convection_coefficient = 7.867410923 # [W/m^2K]
        air_resistance = 1 / (air_convection_coefficient * outer_surface_area)
        return(air_resistance)
    
    def Calculate_Cylinder_Resistance(outer_radius, inner_radius, height, conductivity):
        cylinder_resistance = np.log(outer_radius/inner_radius)/(2*np.pi*height*conductivity)
        return(cylinder_resistance)

    def Calculate_Fluid_Resistance(fluid_convection_coefficient, fluid_inner_surface_area):
        fluid_resistance = 1 / (fluid_convection_coefficient * fluid_inner_surface_area)
        
        return(fluid_resistance)
    
    def Calculate_Total_Resistance(fluid_resistance, tank_wall_resistance, atmosphere_resistance):
        total_resistance = fluid_resistance + tank_wall_resistance + atmosphere_resistance
        
        return(total_resistance)
    
    def Calculate_Rate_of_Heat_Transfer(total_resistance, fluid_temperature, atmosphere_temperature = -20 + 273.15):
        
        # Fluid Temperatures at a Sufficient Distance from the Surface

        rate_of_heat_transfer = (atmosphere_temperature - fluid_temperature)/total_resistance
        
        return(rate_of_heat_transfer)
    
    
    def Calculate_Boil_Off_Rate(rate_of_heat_transfer):
        
        heat_of_fluid_evaporation = 201440.01 # [J/kg]
        boil_off_rate = rate_of_heat_transfer/heat_of_fluid_evaporation
        return(boil_off_rate)
    
    @dataclass
    class Tank:
        tank_name: str
        tank_volume: float
        tank_ID: float
        tank_OD: float
        tank_conductivity: float
        fluid_initial_pressure: float
        fluid_name: str
        fluid_volume: float
        fluid_mass: float = field(init=False)
        
        @property
        def fluid_temperature(self):
            fluid_initial_temperature = PropsSI("T", "P", self.fluid_initial_pressure, "Q", 0, self.fluid_name)
            return(fluid_initial_temperature)
        
        @property
        def fluid_density(self):
            fluid_density = PropsSI("D", "P", self.fluid_initial_pressure, "Q", 0, self.fluid_name)
            return(fluid_density)
        
        @property
        def tank_inner_circular_area(self):
            tank_circular_area = Calculate_Circle_Area(diameter = self.tank_ID)
            return(tank_circular_area)
        
        @property
        def tank_inner_circumference(self):
            tank_inner_circumference =2*np.pi*self.tank_ID
            return(tank_inner_circumference)
        
        @property
        def tank_outer_circumference(self):
            tank_outer_circumference =2*np.pi*self.tank_OD
            return(tank_outer_circumference)
        
        @property
        def fluid_height(self):
            fluid_height = (self.fluid_mass/self.fluid_density) / self.tank_inner_circular_area
            return(fluid_height)
        
        @property
        def fluid_inner_surface_area(self):
            fluid_inner_surface_area = self.tank_inner_circumference * self.fluid_height
            return(fluid_inner_surface_area)
            
        @property
        def tank_outer_surface_area(self):
            tank_outer_surface_area = self.tank_outer_circumference * self.fluid_height
            return(tank_outer_surface_area)
        
        def __post_init__(self):
            self.fluid_mass = self.fluid_density * self.fluid_volume
            
        
    fig, ax = plt.subplots(2, 2, sharex=False)
    
    # for fluid_initial_pressure in [15*c.PSI2PA, 300*c.PSI2PA]:
    for fluid_initial_pressure in [300*c.PSI2PA]:
    
        oxidizer_tank = Tank(
                            tank_name = "oxidizer_tank",
                            tank_volume = parameters.oxidizer_tank_usable_volume,
                            tank_ID = parameters.tank_inner_diameter,
                            tank_OD = parameters.tank_outer_diameter,
                            tank_conductivity = 91.91421797, # aluminum
                            fluid_name = "Oxygen",
                            fluid_initial_pressure = fluid_initial_pressure, #300*c.PSI2PA,
                            fluid_volume = parameters.oxidizer_tank_usable_volume,
                            )
        
        print(f"Tank.fluid_temperature: {oxidizer_tank.fluid_temperature}")
        
        dt = 40
        time = 0

        fluid_mass_array = []
        time_array = []
        rate_of_heat_transfer_array = []
        boil_off_rate_array = []
        fluid_height_array = []
        
        while oxidizer_tank.fluid_mass > 0.05:
            fluid_mass_array.append(oxidizer_tank.fluid_mass)
            time_array.append(time)
            
            time = time + dt
            
            print(f"oxidizer_tank.fluid_mass: {oxidizer_tank.fluid_mass:.9f}")
            
            air_resistance = Calculate_Air_Resistance(outer_surface_area=oxidizer_tank.tank_outer_surface_area)

            tank_wall_resistance = Calculate_Cylinder_Resistance(outer_radius = oxidizer_tank.tank_OD/2, 
                                        inner_radius = oxidizer_tank.tank_ID/2, 
                                        height = oxidizer_tank.fluid_height, 
                                        conductivity = oxidizer_tank.tank_conductivity)

            fluid_resistance = Calculate_Fluid_Resistance(fluid_convection_coefficient = 1000,
                                                        fluid_inner_surface_area = oxidizer_tank.fluid_inner_surface_area) # 10-1000

            total_resistance = Calculate_Total_Resistance(fluid_resistance, 
                                                        tank_wall_resistance, 
                                                        air_resistance)

            rate_of_heat_transfer = Calculate_Rate_of_Heat_Transfer(total_resistance, fluid_temperature = oxidizer_tank.fluid_temperature)
            rate_of_heat_transfer_array.append(rate_of_heat_transfer)
            
            boil_off_rate = Calculate_Boil_Off_Rate(rate_of_heat_transfer)
            boil_off_rate_array.append(boil_off_rate)
            
            oxidizer_tank.fluid_mass -= boil_off_rate * dt
            
            fluid_height_array.append(oxidizer_tank.fluid_height)


        time_array = np.asarray(time_array)
        
        if max(time_array) < 120:
            time_unit_name = "Seconds"
        else:
            time_array = time_array/60
            time_unit_name = "Minutes"    



        ax[0,0].plot(time_array, fluid_mass_array)
        ax[0,0].set_ylabel(f"Liquid Oxygen Mass [kg]")
        ax[0,0].set_xlabel(f"Time [{time_unit_name}]")
        ax[0,0].grid()

        ax[0,1].plot(time_array, boil_off_rate_array)
        ax[0,1].set_ylabel(f"Mass flow rate [kg/s]")
        ax[0,1].set_xlabel(f"Time [{time_unit_name}]")
        ax[0,1].grid()

        # ax[1,0].plot(time_array, pressurant_temperature_array)
        # ax[1,0].set_ylabel(f"Temperature [K]")
        # ax[1,0].set_xlabel(f"Time [{time_unit_name}]")
        # ax[1,0].grid()

        ax[1,1].plot(time_array, fluid_height_array)
        ax[1,1].set_ylabel(f"Fluid Height [m]")
        ax[1,1].set_xlabel(f"Time [{time_unit_name}]")
        ax[1,1].grid()

        # fig.suptitle(f"Emptying {tank_with_orifice.pressurant_name} in {tank_with_orifice.name} through {tank_with_orifice.orifice.nominal_size} orifice using {model_type} model")
   
    plt.show()


def main(parameters):
    parameters = Simulate_Boil_Off(parameters)
    return(parameters)
    

if __name__ == "__main__":
    parameters, wet_mass_distribution, dry_mass_distribution = vehicle_parameters.main()
    main(parameters)
