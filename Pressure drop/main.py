# Pressure drop calculations for Pathfinder Rocket A
# Derived from Rocket 4 code by Keshav Narayanan and Isaiah Jarvis
# Modified and adapted by Luke Goddard

from fluids import fittings
from fluids import core
from CoolProp import CoolProp as CP
from pint import UnitRegistry
import math
import constants
import pandas as pd

u = UnitRegistry()

ox_type = 'oxygen'

# Fluid Properties and other initial parameters
chamber_press = 150 * u.psi  # Pressure in psi
tank_press_ox = 250 * u.psi  # Pressure in psi
tank_press_ipa = 250 * u.psi
mass_flow_ox = 1.78 * (u.pound / u.second)  # LOx mass flow rate
mass_flow_ipa = 1.78 * (u.pound / u.second)  # Isopropyl alcohol mass flow rate

# Conversions into metric units
chamber_press = chamber_press.to(u.pascal)  # Convert to pascals
tank_press_ox = tank_press_ox.to(u.pascal)  # Convert to pascals
tank_press_ipa = tank_press_ipa.to(u.pascal)
mass_flow_ox = mass_flow_ox.to(u.kilogram / u.second)  # Convert to kg/s
mass_flow_ipa = mass_flow_ipa.to(u.kilogram / u.second)  # Convert to kg/s

saturation_temp_ox = CP.PropsSI('T', 'P', tank_press_ox.magnitude, 'Q', 0, ox_type) * u.kelvin
saturation_temp_ipa = 453.01 * u.kelvin # This is sat temp for IPA at 250 psi

density_ox = CP.PropsSI("D", "P", tank_press_ox.magnitude + 10, "T", saturation_temp_ox.magnitude, ox_type) * (u.kilogram / u.meter**3)
dynamic_visc_ox = CP.PropsSI('V', 'P', tank_press_ox.magnitude + 10, 'T', saturation_temp_ox.magnitude, 'oxygen') * (u.pascal * u.second)
kinetic_visc_ox = dynamic_visc_ox / density_ox

density_ipa = constants.DENSITY_IPA * (u.kilogram / u.meter**3)
dynamic_visc_ipa = 0.002055 * (u.pascal * u.second) # No clue if this is right, got from here: https://www.celsius-process.com/wp-content/uploads/2020/01/isopropanol.pdf
kinetic_visc_ipa = dynamic_visc_ipa / density_ipa

# Function by Keshav Narayanan and Isaiah Jarvis from Rocket 4 
def pipe_properties(outer_dia, wall_thic, abs_roughness):
    inner_dia = (outer_dia - 2 * wall_thic).to(u.meter)  # inner diameter

    area = math.pi * (inner_dia / 2)**2   # hydraulic area
    area = area.to(u.meter**2)

    line_vel_ox = (mass_flow_ox / (density_ox * area)).to(u.meter / u.second)    # line velocity
    line_vel_ipa = (mass_flow_ipa / (density_ipa * area)).to(u.meter / u.second)

    rel_rough = core.relative_roughness(inner_dia.magnitude, abs_roughness.to(u.meter).magnitude)   # relative roughness

    reynold_ox = core.Reynolds(D=inner_dia.magnitude, V=line_vel_ox.magnitude, nu=kinetic_visc_ox.magnitude)    # Reynold's number
    reynold_ipa = core.Reynolds(D=inner_dia.magnitude, V=line_vel_ipa.magnitude, nu=kinetic_visc_ipa.magnitude)

    fric_ox = fittings.friction_factor(reynold_ox, eD=rel_rough)    # Darcy friction factor
    fric_ipa = fittings.friction_factor(reynold_ipa, eD=rel_rough)

    return inner_dia, line_vel_ipa, line_vel_ox, reynold_ipa, reynold_ox, fric_ipa, fric_ox, rel_rough

def straight_tube_ox(length, outer_dia, wall_thic, abs_roughness):
    (inner_dia, line_vel_ipa, line_vel_ox, reynold_ipa, reynold_ox, fric_ipa, fric_ox, rel_rough) = pipe_properties(outer_dia, wall_thic, abs_roughness)

    # Calcualation of loss coefficient
    K = core.K_from_f(fric_ox, length.to(u.meter).magnitude, inner_dia.to(u.meter).magnitude)
    drop = core.dP_from_K(K, rho=density_ox.to(u.kilogram / u.meter**3).magnitude, V=line_vel_ox.to(u.meter / u.second).magnitude) * u.pascal
    drop = drop.to(u.psi)
    return drop

def straight_tube_ipa(length, outer_dia, wall_thic, abs_roughness):
    (inner_dia, line_vel_ipa, line_vel_ox, reynold_ipa, reynold_ox, fric_ipa, fric_ox, rel_rough) = pipe_properties(outer_dia, wall_thic, abs_roughness)

    # Calcualation of loss coefficient
    K = core.K_from_f(fric_ipa, length.to(u.meter).magnitude, inner_dia.to(u.meter).magnitude)
    drop = core.dP_from_K(K, rho=density_ipa.to(u.kilogram / u.meter**3).magnitude, V=line_vel_ipa.to(u.meter / u.second).magnitude) * u.pascal
    drop = drop.to(u.psi)
    return drop

def bend_ox(angle, bend_radius, outer_dia, wall_thic, abs_roughness):
    (inner_dia, line_vel_ipa, line_vel_ox, reynold_ipa, reynold_ox, fric_ipa, fric_ox, rel_rough) = pipe_properties(outer_dia, wall_thic, abs_roughness)
    
    bend_radius = bend_radius.to(u.meter)
    # Calculation of loss coefficient
    K = fittings.bend_rounded(
        Di=inner_dia.magnitude,
        angle=angle.magnitude,
        fd=fric_ox,
        rc=bend_radius.magnitude,
        Re=reynold_ox,
        method='Crane'
    )
    drop = core.dP_from_K(K, rho=density_ox.to(u.kilogram / u.meter**3).magnitude, V=line_vel_ox.to(u.meter / u.second).magnitude) * u.pascal
    drop = drop.to(u.psi)
    return drop

def bend_ipa(angle, bend_radius, outer_dia, wall_thic, abs_roughness):
    (inner_dia, line_vel_ipa, line_vel_ox, reynold_ipa, reynold_ox, fric_ipa, fric_ox, rel_rough) = pipe_properties(outer_dia, wall_thic, abs_roughness)
    
    bend_radius = bend_radius.to(u.meter)
    # Calculation of loss coefficient
    K = fittings.bend_rounded(
        Di=inner_dia.magnitude,
        angle=angle.magnitude,
        fd=fric_ipa,
        rc=bend_radius.magnitude,
        Re=reynold_ipa,
        method='Crane'
    )
    drop = core.dP_from_K(K, rho=density_ipa.to(u.kilogram / u.meter**3).magnitude, V=line_vel_ipa.to(u.meter / u.second).magnitude) * u.pascal
    drop = drop.to(u.psi)
    return drop

def ball_valve_ox(cv, valve_dia, wall_thic, abs_roughness):
    (inner_dia, line_vel_ipa, line_vel_ox, reynold_ipa, reynold_ox, fric_ipa, fric_ox, rel_rough) = pipe_properties(valve_dia + (2*wall_thic), wall_thic, abs_roughness)
    
    valve_dia = valve_dia.to(u.meter)

    # Calculation of loss coefficient
    K = fittings.Cv_to_K(Cv=cv, D=valve_dia.magnitude)
    drop = core.dP_from_K(K, rho=density_ox.to(u.kilogram / u.meter**3).magnitude, V=line_vel_ox.to(u.meter / u.second).magnitude) * u.pascal
    drop = drop.to(u.psi)
    return drop

def ball_valve_ipa(cv, valve_dia, wall_thic, abs_roughness):
    (inner_dia, line_vel_ipa, line_vel_ox, reynold_ipa, reynold_ox, fric_ipa, fric_ox, rel_rough) = pipe_properties(valve_dia + (2*wall_thic), wall_thic, abs_roughness)
    
    valve_dia = valve_dia.to(u.meter)

    # Calculation of loss coefficient
    K = fittings.Cv_to_K(Cv=cv, D=valve_dia.magnitude)
    drop = core.dP_from_K(K, rho=density_ipa.to(u.kilogram / u.meter**3).magnitude, V=line_vel_ipa.to(u.meter / u.second).magnitude) * u.pascal
    drop = drop.to(u.psi)
    return drop

def sharp_contraction_ox(outer_dia_one, outer_dia_two, wall_thic, abs_roughness):
    (inner_dia_one, line_vel_ipa, line_vel_ox, reynold_ipa, reynold_ox, fric_ipa, fric_ox, rel_rough) = pipe_properties(outer_dia_one, wall_thic, abs_roughness)
    
    inner_dia_two = (outer_dia_two - 2 * wall_thic).to(u.meter)

    # Calculation of loss coefficient
    K = fittings.contraction_sharp(inner_dia_one, inner_dia_two, fric_ox, reynold_ox, rel_rough, 'Rennels')
    drop = core.dP_from_K(K, rho=density_ox.to(u.kilogram / u.meter**3).magnitude, V=line_vel_ox.to(u.meter / u.second).magnitude) * u.pascal
    drop = drop.to(u.psi)
    return drop

def sharp_contraction_ipa(outer_dia_one, outer_dia_two, wall_thic, abs_roughness):
    (inner_dia_one, line_vel_ipa, line_vel_ox, reynold_ipa, reynold_ox, fric_ipa, fric_ox, rel_rough) = pipe_properties(outer_dia_one, wall_thic, abs_roughness)
    
    inner_dia_two = (outer_dia_two - 2 * wall_thic).to(u.meter)

    # Calculation of loss coefficient
    K = fittings.contraction_sharp(inner_dia_one, inner_dia_two, fric_ipa, reynold_ipa, rel_rough, 'Rennels')
    drop = core.dP_from_K(K, rho=density_ipa.to(u.kilogram / u.meter**3).magnitude, V=line_vel_ipa.to(u.meter / u.second).magnitude) * u.pascal
    drop = drop.to(u.psi)
    return drop

def sharp_expansion_ox(outer_dia_one, outer_dia_two, wall_thic, abs_roughness):
    (inner_dia_one, line_vel_ipa, line_vel_ox, reynold_ipa, reynold_ox, fric_ipa, fric_ox, rel_rough) = pipe_properties(outer_dia_one, wall_thic, abs_roughness)
    
    inner_dia_two = (outer_dia_two - 2 * wall_thic).to(u.meter)

    # Calculation of loss coefficient
    K = fittings.diffuser_sharp(inner_dia_one, inner_dia_two, fric_ox, reynold_ox, rel_rough, 'Rennels')
    drop = core.dP_from_K(K, rho=density_ox.to(u.kilogram / u.meter**3).magnitude, V=line_vel_ox.to(u.meter / u.second).magnitude) * u.pascal
    drop = drop.to(u.psi)
    return drop

def sharp_expansion_ipa(outer_dia_one, outer_dia_two, wall_thic, abs_roughness):
    (inner_dia_one, line_vel_ipa, line_vel_ox, reynold_ipa, reynold_ox, fric_ipa, fric_ox, rel_rough) = pipe_properties(outer_dia_one, wall_thic, abs_roughness)
    
    inner_dia_two = (outer_dia_two - 2 * wall_thic).to(u.meter)

    # Calculation of loss coefficient
    K = fittings.diffuser_sharp(inner_dia_one, inner_dia_two, fric_ipa, reynold_ipa, rel_rough, 'Rennels')
    drop = core.dP_from_K(K, rho=density_ipa.to(u.kilogram / u.meter**3).magnitude, V=line_vel_ipa.to(u.meter / u.second).magnitude) * u.pascal
    drop = drop.to(u.psi)
    return drop
    
def main():
    # Majority of this code is by Keshav Narayanan and Isaiah Jarvis from Rocket 4
    chamber_press = 150 * u.psi  # Pressure in psi
    chamber_press = chamber_press.to(u.pascal)  # Convert to pascals

    injector_inlet_press = chamber_press * 1.2  # Assumes the injector inlet pressure is 20% higher than the chamber pressure

    current_press_ipa = injector_inlet_press    # Set to injector inlet pressure because this is the pressure right before entering injector
    current_press_ox = injector_inlet_press     # Works backwards from the chamber, so this is where we will start. All pressure drops will be added to this number

    reading_ox = False    #The excel sheet will assume it is reading Fuel until the 'Fu' keyword comes up. Then it will switch to Oxygen.

    file = pd.read_excel('pressure_drop_components.xlsx')

    print("\nFuel pressure drops:\n")
    for row in range(len(file)-1, -1, -1):

        name = file.loc[row, 'Name']
        
        if file.loc[row, 'Type'] == 'Fu':
            reading_ox = True
            print("\nLOx pressure drops:\n")

        if file.loc[row, 'Type'] == 'Straight Tube':
            length = file.loc[row, 'Property 1'] * u.inch
            d_outer = file.loc[row, 'Property 3'] * u.inch
            wall_thickness = file.loc[row, 'Property 4'] * u.inch
            abs_rough = file.loc[row, 'Property 5'] * u.mm

            if reading_ox == True:
                part_pd = straight_tube_ox(length, d_outer, wall_thickness, abs_rough)
                print(f"{name}\nPart Type: Straight Tube\nPressure Drop: {part_pd:.2f}\n")
                current_press_ox += part_pd
            else:
                part_pd = straight_tube_ipa(length, d_outer, wall_thickness, abs_rough)
                print(f"{name}\nPart Type: Straight Tube\nPressure Drop: {part_pd:.2f}\n")
                current_press_ipa += part_pd

        if file.loc[row, 'Type'] == 'Bend':
            angle = file.loc[row, 'Property 1'] * u.degree
            bend_radius = file.loc[row, 'Property 2'] * u.inch
            d_outer = file.loc[row, 'Property 3'] * u.inch
            wall_thickness = file.loc[row, 'Property 4'] * u.inch
            abs_rough = file.loc[row, 'Property 5'] * u.mm

            if reading_ox == True:
                part_pd = bend_ox(angle, bend_radius, d_outer, wall_thickness, abs_rough)
                print(f"{name}\nPart Type: Bend\nPressure Drop: {part_pd:.2f}\n")
                current_press_ox += part_pd
            else:
                part_pd = bend_ipa(angle, bend_radius, d_outer, wall_thickness, abs_rough)
                print(f"{name}\nPart Type: Bend\nPressure Drop: {part_pd:.2f}\n")
                current_press_ipa += part_pd


        if file.loc[row, 'Type'] == 'Ball Valve':
            inner_dia = file.loc[row, 'Property 1'] * u.inch
            Cv = file.loc[row, 'Property 2']
            d_outer = file.loc[row, 'Property 3'] * u.inch
            wall_thickness = file.loc[row, 'Property 4'] * u.inch
            abs_rough = file.loc[row, 'Property 5'] * u.mm

            if reading_ox == True:
                part_pd = ball_valve_ox(Cv, inner_dia, wall_thickness, abs_rough)
                print(f"{name}\nPart Type: Ball Valve\nPressure Drop: {part_pd:.2f}\n")
                current_press_ox += part_pd
            else:
                part_pd = ball_valve_ipa(Cv, inner_dia, wall_thickness, abs_rough)
                print(f"{name}\nPart Type: Ball Valve\nPressure Drop: {part_pd:.2f}\n")
                current_press_ipa += part_pd

        if file.loc[row, 'Type'] == 'Venturi':
            if reading_ox == True:
                venturi_inlet_press = current_press_ox / .8    #assumes an 80% pressure recovery rate
                venturi_pd = venturi_inlet_press - current_press_ox
                print(f"{name}\nPart Type: Venturi\nPressure Drop: {venturi_pd:.2f}\n")
                current_press_ox += venturi_pd
            else:
                venturi_inlet_press = current_press_ipa / .8
                venturi_pd = venturi_inlet_press - current_press_ipa
                print(f"{name}\nPart Type: Venturi\nPressure Drop: {venturi_pd:.2f}\n")
                current_press_ipa += venturi_pd

            


    print("\nTotal system pressure drops:")
    print(f"Oxygen system: {(current_press_ox - chamber_press).to(u.psi):.2f}")
    print(f"Ethanol system: {(current_press_ipa - chamber_press).to(u.psi):.2f}")
    print(f"With a LOx chamber pressure of {chamber_press.to(u.psi):.2f}, the tank pressure will be {current_press_ox.to(u.psi):.2f}")
    print(f"With an ethanol chamber pressure of {chamber_press.to(u.psi):.2f}, the tank pressure will be {current_press_ipa.to(u.psi):.2f}")

if __name__ == "__main__":
    main()