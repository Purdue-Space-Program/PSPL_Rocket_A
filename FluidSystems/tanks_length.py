import os
import sys
import matplotlib.pyplot as plt
import numpy as np
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
from vehicle_parameters import parameters as p

def calc_ullage_volume_bulkhead(radius, depth):
    return ((np.pi * depth**2 / 3) * (3 * radius - depth))

def calc_residual_fraction_from_AoA(diameter, angle, fuel_use_volume):
    height = np.tan(np.radians(angle)) * diameter
    volume = 0.5 * np.pi * diameter**2 / 4 * height
    fraction = volume / fuel_use_volume
    return fraction

##### FUEL TANK PARAMETERS #####
fuel_tank_ID = p.tube_inner_diameter
fuel_use_volume = p.fuel_tank_usable_volume
fuel_ullage_fraction = 10 / 100

angle = 0
angle_history = []
fraction_history = []
while angle < 45:
    residual_fraction = calc_residual_fraction_from_AoA(fuel_tank_ID, angle, fuel_use_volume)
    fraction_history.append(residual_fraction)
    angle_history.append(angle)
    angle +=1

plt.plot(angle_history, np.array(fraction_history) * 100)
plt.xlabel("AoA [Degrees]")
plt.ylabel("Residual %")
plt.grid()
plt.title("AoA vs Residual Fraction")
plt.show()


fuel_residual_fraction = 20 / 100

##### LOX TANK PARAMETERS #####
ox_tank_ID = p.tube_inner_diameter
ox_use_volume = p.oxidizer_tank_usable_volume
ox_ullage_fraction = 10 / 100
ox_residual_fraction = 20 / 100

##### BULKHEAD PARAMETERS #####
bulkhead_radius_of_curvature = 3.5 * c.IN2M
bulkhead_depth = 1 * c.IN2M
bulkhead_tank_slot = 1.1 * c.IN2M # Length of tank overlapping with outer surface of bulkhead 

##### Calculate total volume #####
bulkhead_ullage_volume = calc_ullage_volume_bulkhead(bulkhead_radius_of_curvature, bulkhead_depth)

fuel_ullage_volume = fuel_use_volume * fuel_ullage_fraction + bulkhead_ullage_volume
fuel_residual_volume = fuel_use_volume * fuel_residual_fraction
fuel_tank_volume = fuel_ullage_volume + fuel_residual_volume + fuel_use_volume

ox_ullage_volume = ox_use_volume * ox_ullage_fraction + bulkhead_ullage_volume
ox_residual_volume = ox_use_volume * ox_residual_fraction
ox_tank_volume = ox_ullage_volume + ox_residual_volume + ox_use_volume

##### Calculate total length #####
fuel_tank_area = np.pi / 4 * fuel_tank_ID**2
fuel_tank_length = (((fuel_ullage_volume - bulkhead_ullage_volume) + fuel_use_volume + fuel_residual_volume) / fuel_tank_area) + bulkhead_tank_slot * 2

ox_tank_area = np.pi / 4 * ox_tank_ID**2
ox_tank_length = (((ox_ullage_volume - bulkhead_ullage_volume) + ox_use_volume + ox_residual_volume) / ox_tank_area) + bulkhead_tank_slot * 2


print("-----PARAMETERS-----")
print(f"Fuel Ullage: {fuel_ullage_fraction * 100} %")
print(f"Fuel Residual: {fuel_residual_fraction * 100} %")
print(f"Ox Ullage: {ox_ullage_fraction * 100} %")
print(f"Ox Residual: {ox_residual_fraction * 100} %")
print("-----RESULTS-----")
print(f"Fuel Tube Length: {fuel_tank_length * c.M2IN:.1f} in")
print(f"Ox Tube Length: {ox_tank_length * c.M2IN:.1f} in")
print(f"Fuel Tank Length: {(fuel_tank_length - bulkhead_tank_slot * 2) * c.M2IN:.1f} in")
print(f"Ox Tank Length {(ox_tank_length - bulkhead_tank_slot * 2) * c.M2IN:.1f} in")
print(f"Real Fuel Ullage: {fuel_ullage_volume / fuel_tank_volume * 100:.1f} %")
print(f"Real Ox Ullage: {ox_ullage_volume / ox_tank_volume * 100:.1f} %")