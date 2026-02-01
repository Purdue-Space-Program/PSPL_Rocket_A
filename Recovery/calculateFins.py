import sys
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *
from vehicle_parameters import parameters
import constants as c
import csv

def calculateCenterOfPressure(target_cal, center_of_mass, diameter):
    return ((target_cal)*(diameter))+center_of_mass

def calculateTaperRatio(tip_chord, root_chord):
    return tip_chord/root_chord

def calculateCriticalMachNumber(shear_modulus, aspect_ratio, fin_thickness, c, C, pressure, taper_ratio,):
    num = 2 * shear_modulus * (aspect_ratio + 2) * (fin_thickness/fin_thickness) * (fin_thickness/fin_thickness) * (fin_thickness/fin_thickness)
    dem = C * aspect_ratio * aspect_ratio * aspect_ratio * pressure * (taper_ratio + 1)
    return np.sqrt(num/dem) 

def main():
    target_stability_caliber = 2.0
    wet_target_center_of_pressure_from_top = calculateCenterOfPressure(target_stability_caliber, parameters.wet_COM_location_from_top, parameters.tube_outer_diameter)
    dry_target_center_of_pressure_from_top = calculateCenterOfPressure(target_stability_caliber, parameters.dry_COM_location_from_top, parameters.tube_outer_diameter)

    print(f"wet_target_center_of_pressure_from_top: {wet_target_center_of_pressure_from_top * c.M2IN:.2f} inches from top")
    print(f"dry_target_center_of_pressure_from_top: {dry_target_center_of_pressure_from_top * c.M2IN:.2f} inches from top")

main()