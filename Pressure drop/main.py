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

def sharp_edged_inlet():
    K = 0.5
    return K

def inward_projecting_sharp_edged_inlet():
    K = 0.78
    return K

def sharp_edged_inlet_at_angle(angle):
    K = 0.5 + 0.3 * math.cos(angle) + 0.2 * math.cos(angle)**2
    return K

def rounded_inlet(radius, diameter):
    r_by_d = radius / diameter
    if r_by_d > 0.2:
        K = 0.03
    else:
        K = 0.5 - 1.07 * r_by_d**(0.5) - 2.13 * r_by_d + 8.24 * r_by_d**(1.5) - 8.48 * r_by_d**2 + 2.90 * r_by_d**(2.5)
    return K

def sharp_edged_sudden_contraction(diameter_one, diameter_two):
    K = 0.5 * (1 - (diameter_one / diameter_two)**2)
    return K

def sharp_edged_sudden_expansion(diameter_one, diameter_two):
    K = 1.05 * (1 - (diameter_one / diameter_two)**2)**2
    return K

def outlet():
    K = 1.05
    return K

def bend(angle, radius, diameter, friction_factor):
    r_by_d = radius / diameter
    K = friction_factor * angle * r_by_d + (0.1 + 2.4 * friction_factor) * math.sin(angle / 2) + (6.6 * friction_factor * math.sin(angle / 2)**0.5 + math.sin(angle / 2)) / (r_by_d**((4 * angle) / math.pi))
    return K

def valve(flow_coefficient, diameter):
    K = (890.4 * diameter**4) / (flow_coefficient**2)
    return K

def main():
    # Majority of this code is by Keshav Narayanan and Isaiah Jarvis from Rocket 4
    reading_ox = False    #The excel sheet will assume it is reading Fuel until the 'Fu' keyword comes up. Then it will switch to Oxygen.

    file = pd.read_excel('pressure_drop_components.xlsx')

    print("\nFuel pressure drops:\n")
    for row in range(len(file)-1, -1, -1):

        name = file.loc[row, 'Name']
        
        if file.loc[row, 'Type'] == 'Fu':
            reading_ox = True
            print("\nLOx pressure drops:\n")

        if file.loc[row, 'Type'] == 'Straight Tube':

            '''
            length = file.loc[row, 'Property 1'] * u.inch
            d_outer = file.loc[row, 'Property 3'] * u.inch
            wall_thickness = file.loc[row, 'Property 4'] * u.inch
            abs_rough = file.loc[row, 'Property 5'] * u.mm
            '''

            if reading_ox == True:
                #part_pd = straight_dp_ox(length, d_outer, wall_thickness, abs_rough)
                print(f"{name}\nPart Type: Straight Tube\nPressure Drop: {1:.2f}\n")
                #current_press_ox += part_pd
            else:
                #part_pd = straight_dp_eth(length, d_outer, wall_thickness, abs_rough)
                print(f"{name}\nPart Type: Straight Tube\nPressure Drop: {1:.2f}\n")
                #current_press_eth += part_pd

        if file.loc[row, 'Type'] == 'Bend':

            '''angle = file.loc[row, 'Property 1'] * u.degree
            bend_radius = file.loc[row, 'Property 2'] * u.inch
            d_outer = file.loc[row, 'Property 3'] * u.inch
            wall_thickness = file.loc[row, 'Property 4'] * u.inch
            abs_rough = file.loc[row, 'Property 5'] * u.mm'''

            if reading_ox == True:
                #part_pd = bend_dp_ox(angle, bend_radius, d_outer, wall_thickness, abs_rough)
                print(f"{name}\nPart Type: Bend\nPressure Drop: {1:.2f}\n")
                #current_press_ox += part_pd
            else:
                #part_pd = bend_dp_eth(angle, bend_radius, d_outer, wall_thickness, abs_rough)
                print(f"{name}\nPart Type: Bend\nPressure Drop: {1:.2f}\n")
                #current_press_eth += part_pd


        if file.loc[row, 'Type'] == 'Ball Valve':

            '''inner_dia = file.loc[row, 'Property 1'] * u.inch
            Cv = file.loc[row, 'Property 2']
            d_outer = file.loc[row, 'Property 3'] * u.inch
            wall_thickness = file.loc[row, 'Property 4'] * u.inch
            abs_rough = file.loc[row, 'Property 5'] * u.mm'''

            if reading_ox == True:
                #part_pd = ball_valve_dp_ox(Cv, inner_dia, wall_thickness, abs_rough)
                print(f"{name}\nPart Type: Ball Valve\nPressure Drop: {1:.2f}\n")
                #current_press_ox += part_pd
            else:
                #part_pd = ball_valve_dp_eth(Cv, inner_dia, wall_thickness, abs_rough)
                print(f"{name}\nPart Type: Ball Valve\nPressure Drop: {1:.2f}\n")
                #current_press_eth += part_pd

        if file.loc[row, 'Type'] == 'Venturi':
            if reading_ox == True:
                #venturi_inlet_press = current_press_ox / .8    #assumes an 80% pressure recovery rate
                #venturi_pd = venturi_inlet_press - current_press_ox
                print(f"{name}\nPart Type: Venturi\nPressure Drop: {1:.2f}\n")
                #current_press_ox += venturi_pd
            else:
                #venturi_inlet_press = current_press_eth / .8
                #venturi_pd = venturi_inlet_press - current_press_eth
                print(f"{name}\nPart Type: Venturi\nPressure Drop: {1:.2f}\n")
                #current_press_eth += venturi_pd

            


    print("\nTotal system pressure drops:")
    '''
    print(f"Oxygen system: {(current_press_ox - chamber_press).to(u.psi):.2f}")
    print(f"Ethanol system: {(current_press_eth - chamber_press).to(u.psi):.2f}")
    print(f"With a LOx chamber pressure of {chamber_press.to(u.psi):.2f}, the tank pressure will be {current_press_ox.to(u.psi):.2f}")
    print(f"With an ethanol chamber pressure of {chamber_press.to(u.psi):.2f}, the tank pressure will be {current_press_eth.to(u.psi):.2f}")
    '''

if __name__ == "__main__":
    main()