# Pressure drop calculations for Pathfinder Rocket A
# Derived from Rocket 4 code by Keshav Narayanan and Isaiah Jarvis
# Modified and adapted by Luke Goddard

from fluids import fittings
from fluids import core
from CoolProp import CoolProp as CP
from pint import UnitRegistry
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import constants as c
from vehicle_parameters import parameters

u = UnitRegistry()

oxidizer_name = "oxygen"

# Fluid Properties and other initial parameters
chamber_pressure =  parameters.chamber_pressure * u.pascal
oxidizer_tank_pressure = parameters.tank_pressure * u.pascal
fuel_tank_pressure = parameters.tank_pressure * u.pascal
fuel_mass_flow_rate = parameters.core_fuel_mass_flow_rate * (u.kilogram / u.second)  # Isopropyl alcohol mass flow rate
oxidizer_mass_flow_rate = parameters.oxidizer_mass_flow_rate * (u.kilogram / u.second)  # LOx mass flow rate

saturation_temp_ox = CP.PropsSI('T', 'P', oxidizer_tank_pressure.magnitude, 'Q', 0, oxidizer_name) * u.kelvin
saturation_temp_ipa = 453.01 * u.kelvin # This is sat temp for IPA at 250 psi

oxidizer_density = CP.PropsSI("D", "P", oxidizer_tank_pressure.magnitude + 10, "T", saturation_temp_ox.magnitude, oxidizer_name) * (u.kilogram / u.meter**3)
oxidizer_dynamic_viscosity = CP.PropsSI('V', 'P', oxidizer_tank_pressure.magnitude + 10, 'T', saturation_temp_ox.magnitude, oxidizer_name) * (u.pascal * u.second)
oxidizer_kinetic_viscosity = oxidizer_dynamic_viscosity / oxidizer_density

fuel_density = c.DENSITY_IPA * (u.kilogram / u.meter**3)
fuel_dynamic_viscosity = 0.002055 * (u.pascal * u.second) # No clue if this is right, got from here: https://www.celsius-process.com/wp-content/uploads/2020/01/isopropanol.pdf
fuel_kinetic_viscosity = fuel_dynamic_viscosity / fuel_density

venturi_pressure_drop_percent = 0.8
injector_pressure_drop_percent = 0.8


def CalculatePressureDropFromFrictionFactor(fluid_friction_factor, tube_length, tube_inner_diameter, fluid_density, fluid_line_velocity):
    K = core.K_from_f(fluid_friction_factor, tube_length.to(u.meter).magnitude, tube_inner_diameter.to(u.meter).magnitude)
    pressure_drop = core.dP_from_K(K, rho=fluid_density.to(u.kilogram / u.meter**3).magnitude, V=fluid_line_velocity.to(u.meter / u.second).magnitude) * u.pascal
    pressure_drop = pressure_drop.to(u.psi)
    return (pressure_drop)

# Function by Keshav Narayanan and Isaiah Jarvis from Rocket 4 
def CalculateFluidInPipe(tube_outer_diameter, tube_wall_thickness, absolute_roughness, fluid_mass_flow_rate, fluid_kinetic_viscosity, fluid_density):
    tube_inner_diameter = (tube_outer_diameter - 2 * tube_wall_thickness).to(u.meter)  # inner diameter

    area = np.pi * (tube_inner_diameter / 2)**2 # hydraulic area
    area = area.to(u.meter**2)

    fluid_line_velocity = (fluid_mass_flow_rate / (fluid_density * area)).to(u.meter / u.second)
    rel_rough = core.relative_roughness(tube_inner_diameter.magnitude, absolute_roughness.to(u.meter).magnitude)   # relative roughness

    fluid_Reynolds_number = core.Reynolds(D=tube_inner_diameter.magnitude, V=fluid_line_velocity.magnitude, nu=fluid_kinetic_viscosity.magnitude)    # Reynold's number

    fluid_friction_factor = fittings.friction_factor(fluid_Reynolds_number, eD=rel_rough)    # Darcy friction factor

    return (tube_inner_diameter, fluid_line_velocity, fluid_Reynolds_number, fluid_friction_factor, rel_rough)

def CalculateStraightTubePressureDrop(tube_length, tube_outer_diameter, tube_wall_thickness, absolute_roughness, fluid_mass_flow_rate, fluid_kinetic_viscosity, fluid_density):
    (tube_inner_diameter, 
     fluid_line_velocity, 
     fluid_Reynolds_number, 
     fluid_friction_factor, 
     rel_rough) = CalculateFluidInPipe(tube_outer_diameter, tube_wall_thickness, absolute_roughness, fluid_mass_flow_rate, fluid_kinetic_viscosity, fluid_density)

    pressure_drop = CalculatePressureDropFromFrictionFactor(fluid_friction_factor, tube_length, tube_inner_diameter, fluid_density, fluid_line_velocity)
    return pressure_drop


def CalculateTubeBendPressureDrop(angle, bend_radius, tube_outer_diameter, tube_wall_thickness, absolute_roughness, fluid_mass_flow_rate, fluid_kinetic_viscosity, fluid_density):
    (tube_inner_diameter, fluid_line_velocity, fluid_Reynolds_number, fluid_friction_factor, rel_rough) = CalculateFluidInPipe(tube_outer_diameter, tube_wall_thickness, absolute_roughness, fluid_mass_flow_rate, fluid_kinetic_viscosity, fluid_density)
    
    bend_radius = bend_radius.to(u.meter)
    # Calculation of loss coefficient
    K = fittings.bend_rounded(
        Di=tube_inner_diameter.magnitude,
        angle=angle.magnitude,
        fd=fluid_friction_factor,
        rc=bend_radius.magnitude,
        bend_diameters = None,
        Re=fluid_Reynolds_number,
        method='Crane'
    )
    

    drop = core.dP_from_K(K, rho=fluid_density.to(u.kilogram / u.meter**3).magnitude, V=fluid_line_velocity.to(u.meter / u.second).magnitude) * u.pascal
    drop = drop.to(u.psi)
    return drop

def CalculateBallValvePressureDrop(Cv, valve_inner_diameter, tube_wall_thickness, absolute_roughness, fluid_mass_flow_rate, fluid_kinetic_viscosity, fluid_density):
    tube_outer_diameter = valve_inner_diameter + (2*tube_wall_thickness)
    (inner_dia, fluid_line_velocity, fluid_Reynolds_number, fluid_friction_factor, rel_rough) = CalculateFluidInPipe(tube_outer_diameter, tube_wall_thickness, absolute_roughness, fluid_mass_flow_rate, fluid_kinetic_viscosity, fluid_density)
        
    valve_inner_diameter = valve_inner_diameter.to(u.meter)

    K = fittings.Cv_to_K(Cv=Cv, D=valve_inner_diameter.magnitude)
    pressure_drop = core.dP_from_K(K, rho=fluid_density.to(u.kilogram / u.meter**3).magnitude, V=fluid_line_velocity.to(u.meter / u.second).magnitude) * u.pascal
    pressure_drop = pressure_drop.to(u.psi)
    return pressure_drop

# def sharp_contraction_ox(outer_dia_one, outer_dia_two, wall_thickness, abs_roughness):
#     (inner_dia_one, fuel_line_velocity, oxidizer_line_velocity, reynold_ipa, reynold_ox, friction_factor_ipa, friction_factor_oxidizer, rel_rough) = CalculateFluidInPipe(outer_dia_one, wall_thickness, abs_roughness)
    
#     inner_dia_two = (outer_dia_two - 2 * wall_thickness).to(u.meter)

#     # Calculation of loss coefficient
#     K = fittings.contraction_sharp(inner_dia_one, inner_dia_two, friction_factor_oxidizer, reynold_ox, rel_rough, "Rennels")
#     drop = core.dP_from_K(K, rho=oxidizer_density.to(u.kilogram / u.meter**3).magnitude, V=oxidizer_line_velocity.to(u.meter / u.second).magnitude) * u.pascal
#     drop = drop.to(u.psi)
#     return drop

# def sharp_contraction_ipa(outer_dia_one, outer_dia_two, wall_thickness, abs_roughness):
#     (inner_diameter_one, fuel_line_velocity, oxidizer_line_velocity, reynold_ipa, reynold_ox, friction_factor_ipa, friction_factor_oxidizer, rel_rough) = CalculateFluidInPipe(outer_dia_one, wall_thickness, abs_roughness)
    
#     inner_diameter_two = (outer_dia_two - 2 * wall_thickness).to(u.meter)

#     # Calculation of loss coefficient
#     K = fittings.contraction_sharp(inner_diameter_one, inner_diameter_two, friction_factor_ipa, reynold_ipa, rel_rough, "Rennels")
#     drop = core.dP_from_K(K, rho=fuel_density.to(u.kilogram / u.meter**3).magnitude, V=fuel_line_velocity.to(u.meter / u.second).magnitude) * u.pascal
#     drop = drop.to(u.psi)
#     return drop

# def sharp_expansion_ox(outer_dia_one, outer_dia_two, wall_thickness, abs_roughness):
#     (inner_dia_one, fuel_line_velocity, oxidizer_line_velocity, reynold_ipa, reynold_ox, friction_factor_ipa, friction_factor_oxidizer, rel_rough) = CalculateFluidInPipe(outer_dia_one, wall_thickness, abs_roughness)
    
#     inner_dia_two = (outer_dia_two - 2 * wall_thickness).to(u.meter)

#     # Calculation of loss coefficient
#     K = fittings.diffuser_sharp(inner_dia_one, inner_dia_two, friction_factor_oxidizer, reynold_ox, rel_rough, "Rennels")
#     drop = core.dP_from_K(K, rho=oxidizer_density.to(u.kilogram / u.meter**3).magnitude, V=oxidizer_line_velocity.to(u.meter / u.second).magnitude) * u.pascal
#     drop = drop.to(u.psi)
#     return drop

# def sharp_expansion_ipa(outer_dia_one, outer_dia_two, wall_thickness, abs_roughness):
#     (inner_dia_one, fuel_line_velocity, oxidizer_line_velocity, reynold_ipa, reynold_ox, friction_factor_ipa, friction_factor_oxidizer, rel_rough) = CalculateFluidInPipe(outer_dia_one, wall_thickness, abs_roughness)
    
#     inner_dia_two = (outer_dia_two - 2 * wall_thickness).to(u.meter)

#     # Calculation of loss coefficient
#     K = fittings.diffuser_sharp(inner_dia_one, inner_dia_two, friction_factor_ipa, reynold_ipa, rel_rough, "Rennels")
#     drop = core.dP_from_K(K, rho=fuel_density.to(u.kilogram / u.meter**3).magnitude, V=fuel_line_velocity.to(u.meter / u.second).magnitude) * u.pascal
#     drop = drop.to(u.psi)
#     return drop
    

# Majority of this code is by Keshav Narayanan and Isaiah Jarvis from Rocket 4
try:
    file = pd.read_excel("Pressure drop/fluid_components.xlsx")
except:
    file = pd.read_excel("fluid_components.xlsx")

number_of_parts = len(file)

# Works backwards from the chamber, so this is where we will start. All pressure drops will be added to this number
current_fuel_pressure = chamber_pressure
current_oxidizer_pressure = chamber_pressure

fuel_pressure_array = [current_fuel_pressure.to(u.psi)]
oxidizer_pressure_array = [current_oxidizer_pressure.to(u.psi)]

fuel_part_names = ["Fuel in Chamber"]
oxidizer_part_names = ["Oxidizer in Chamber"]


reading_ox = True    #The excel sheet will assume it is reading Ox until the 'Fuel' keyword comes up. Then it will switch to Fuel.
print("\n--------LOx pressure drops--------:")

for row in range(number_of_parts-1, -1, -1):

    part_name = file.loc[row, "Part Name"]
    part_type = file.loc[row, "Part Type"]
    
    is_nan = isinstance(part_name, (int, float, np.floating)) and np.isnan(part_name)
    
    if part_name == "Fuel":
        # permanently switch to fuel
        reading_ox = False
        print("\n--------Fuel pressure drops--------:")

    elif isinstance(part_name, str) and part_name.strip() != "":
        tube_length = file.loc[row, "Tube Length [in]"] * u.inch
        inner_diameter = file.loc[row, "Inner Diameter [in]"] * u.inch
        tube_outer_diameter = file.loc[row, "Outer Diameter [in]"] * u.inch
        tube_wall_thickness = file.loc[row, "Wall Thickness [in]"] * u.inch
        absolute_roughness = file.loc[row, 'Absolute Roughness [mm]'] * u.mm
        angle = file.loc[row, "Angle [degrees]"] * u.degree
        bend_radius = file.loc[row, "Bend Radius [in]"] * u.inch
        Cv = file.loc[row, 'Cv']
        

        if reading_ox == True:
            fluid_mass_flow_rate = oxidizer_mass_flow_rate
            fluid_kinetic_viscosity = oxidizer_kinetic_viscosity
            fluid_density = oxidizer_density
            current_fluid_pressure = current_oxidizer_pressure
        else:
            fluid_mass_flow_rate = fuel_mass_flow_rate
            fluid_kinetic_viscosity = fuel_kinetic_viscosity
            fluid_density = fuel_density
            current_fluid_pressure = current_fuel_pressure

        match part_type:
            case "Straight Tube":
                part_pressure_drop = CalculateStraightTubePressureDrop(tube_length, tube_outer_diameter, tube_wall_thickness, absolute_roughness, fluid_mass_flow_rate, fluid_kinetic_viscosity, fluid_density)
            
            case "Bend":
                part_pressure_drop = CalculateTubeBendPressureDrop(angle, bend_radius, tube_outer_diameter, tube_wall_thickness, absolute_roughness, fluid_mass_flow_rate, fluid_kinetic_viscosity, fluid_density)

            case "Ball Valve":
                part_pressure_drop = CalculateBallValvePressureDrop(Cv, inner_diameter, tube_wall_thickness, absolute_roughness, fluid_mass_flow_rate, fluid_kinetic_viscosity, fluid_density)    

            case "Venturi":
                venturi_inlet_pressure = current_fluid_pressure / venturi_pressure_drop_percent
                part_pressure_drop = venturi_inlet_pressure - current_fluid_pressure
            
            case "Injector":
                injector_inlet_pressure = current_fluid_pressure / injector_pressure_drop_percent
                part_pressure_drop = injector_inlet_pressure - current_fluid_pressure

            case _:
                raise ValueError(f"Unknown part type: {part_name}")

        print(f"Part Type: {part_type}\nPart Name: {part_name}\nPressure Drop: {part_pressure_drop.to(u.psi):.2f}\n")



        if reading_ox == True:
            current_oxidizer_pressure += part_pressure_drop
            oxidizer_part_names.append(part_name)
            oxidizer_pressure_array.append(current_oxidizer_pressure.to(u.psi))
        else:
            current_fuel_pressure += part_pressure_drop
            fuel_pressure_array.append(current_fuel_pressure.to(u.psi))
            fuel_part_names.append(part_name)
        
            
    elif is_nan:
        print("NaN!")
    
    else:
        raise ValueError("help")


oxidizer_total_pressure_drop = current_oxidizer_pressure - chamber_pressure
fuel_total_pressure_drop = current_fuel_pressure - chamber_pressure


print(f"IPA: \n\tTotal pressure drop: {fuel_total_pressure_drop.to(u.psi):.2f}\n\tChamber pressure: {chamber_pressure.to(u.psi):.2f}\n\tTank pressure: {current_fuel_pressure.to(u.psi):.2f}")
print(f"LOx: \n\tTotal pressure drop: {oxidizer_total_pressure_drop.to(u.psi):.2f}\n\tChamber pressure: {chamber_pressure.to(u.psi):.2f}\n\tTank pressure: {current_oxidizer_pressure.to(u.psi):.2f}")

if __name__ == "__main__":
    
    # go from left to right
    fuel_part_names.reverse()
    fuel_pressure_array.reverse()
    oxidizer_part_names.reverse()
    oxidizer_pressure_array.reverse()
    
    
    
    use_subplots = True
    
    for count, fluid_pressure_array in enumerate([fuel_pressure_array, oxidizer_pressure_array]):
        
        if fluid_pressure_array is fuel_pressure_array:
            fluid_part_names = fuel_part_names
            fluid_color = "r"
            fluid_name = "Fuel"
        elif fluid_pressure_array is oxidizer_pressure_array:
            fluid_part_names = oxidizer_part_names
            fluid_color = "b"
            fluid_name = "Oxidizer"
        
        
        if use_subplots == True:
            plt.subplot(1, 2, count+1)
            
            # make go in between ticks
            fluid_point_indices = np.arange(len(fluid_pressure_array))          
            fluid_label_positions = fluid_point_indices[:-1] + 0.5
            plt.xlim(-0.5, len(fluid_point_indices) - 0.5)
            plt.xticks(
                ticks=fluid_label_positions,
                labels=fluid_part_names[:-1],
                rotation=45,
                ha="right"
            )  
            plt.grid(True)
            plt.title(fluid_name)
        else:
            plt.legend(["Fuel", "Oxidizer"])


        
        plt.plot(fluid_point_indices, list(float(fluid_pressure.to(u.psi).magnitude) for fluid_pressure in fluid_pressure_array), fluid_color)
        plt.xlabel("Part Name")
        plt.ylabel("Fluid Pressure [psi]")
    
    plt.show()