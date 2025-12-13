'''

Welcome to GDB Online.
GDB online is an online compiler and debugger tool for C, C++, Python, Java, PHP, Ruby, Perl,
C#, OCaml, VB, Swift, Pascal, Fortran, Haskell, Objective-C, Assembly, HTML, CSS, JS, SQLite, Prolog.
Code, Compile, Run and Debug online from anywhere in world.

'''

import numpy as np

import math
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
import vehicle_parameters as vehicle

def main():
    max_height = vehicle.parameters.estimated_apogee # [m]
    rocket_mass = vehicle.parameters.dry_mass # [kg]
    air_density = 1.229 # [kg/m^3] https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/rktvrecv.html

    parachute_name = "Fruity Chutes"
    # parachute_name = "Spherachutes"

    if parachute_name == "Fruity Chutes":
        drag_coefficient = 2.2 # https://shop.fruitychutes.com/collections/parachutes/products/iris-ultra-120-standard-parachute-79lbs-20fps
        # effective_parachute_area = 24.6677824063 abhi?
        outer_parachute_diameter = 192 * c.IN2M # https://shop.fruitychutes.com/collections/parachutes/products/iris-ultra-120-standard-parachute-79lbs-20fps
        spill_hole_diameter = 15.76 * c.IN2M # guess
        effective_parachute_area = CalculateEffectiveParachuteArea(outer_parachute_diameter, spill_hole_diameter)
          
        outer_parachute_diameter = 122.23 * c.IN2M # https://spherachutes.com/products/192-inch-spherachute
        spill_hole_diameter = 15.76 * c.IN2M # https://spherachutes.com/products/192-inch-spherachute
    elif parachute_name == "Spherachutes":
        drag_coefficient = 0.75 # https://spherachutes.com/pages/decent-rate-chart
        outer_parachute_diameter = 122.23 * c.IN2M # https://spherachutes.com/products/192-inch-spherachute
        spill_hole_diameter = 15.76 * c.IN2M # https://spherachutes.com/products/192-inch-spherachute
        effective_parachute_area = CalculateEffectiveParachuteArea(outer_parachute_diameter, spill_hole_diameter)
        print(f"effective_parachute_area: {effective_parachute_area:.2f}")
        # canopy_area = 15.029593683 abhi?
    else:
        raise ValueError("!!parachute does not exist!!")

    weight_force = c.GRAVITY * rocket_mass
    terminal_velocity = math.sqrt((2*weight_force)/(effective_parachute_area*drag_coefficient*air_density)) # solving for velocity setting weight and Drag equal
    drag_force = 0.5 * drag_coefficient * air_density * effective_parachute_area *(terminal_velocity**2) # https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/rktvrecv.html
    descent_time = max_height/terminal_velocity

    print("\nVehicle:")
    print(f"\tMass:{rocket_mass * c.KG2LBM:.2f} lbm")
    print(f"\tDrag Force: {drag_force:.2f} N")
    print(f"\tWeight Force: {weight_force:.2f} N")
    
    print(f"\nParachute: {parachute_name}")
    print(f"\tCd: {drag_coefficient:.2f}")
    print(f"\tArea: {effective_parachute_area:.2f} m^2")
    print(f"\tTerminal Velocity: {terminal_velocity * c.M2FT:.2f} ft/s")
    print(f"\tDescent Time: {descent_time:.2f} seconds")

# account for parachute spill hole
def CalculateEffectiveParachuteArea(outer_diameter, spill_hole_diameter):
    total_area = CalculateCircleArea(outer_diameter)
    spill_hole_area = CalculateCircleArea(spill_hole_diameter)
    effective_parachute_area = total_area - spill_hole_area
    return(effective_parachute_area)
    
def CalculateCircleArea(diameter):
    radius = diameter/2
    area = np.pi * (radius**2)
    return(area)

if __name__ == "__main__":
    main()