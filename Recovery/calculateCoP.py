import sys
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *
from vehicle_parameters import parameters
import constants as c
import csv

def calculateCoP(target_cal, CoM, diameter):
    return ((target_cal)*(diameter))+CoM

def main():
    target_stability_caliber = 2.0
    wet_target_center_of_pressure_from_top = calculateCoP(target_stability_caliber, parameters.wet_COM_location_from_top, parameters.tube_outer_diameter)
    dry_target_center_of_pressure_from_top = calculateCoP(target_stability_caliber, parameters.dry_COM_location_from_top, parameters.tube_outer_diameter)

    print(f"wet_target_center_of_pressure_from_top: {wet_target_center_of_pressure_from_top * c.M2IN:.2f} inches from top")
    print(f"dry_target_center_of_pressure_from_top: {dry_target_center_of_pressure_from_top * c.M2IN:.2f} inches from top")

main()