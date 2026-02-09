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

def calculateStabilityCaliber(center_of_pressure, center_of_mass, diameter):
    return (center_of_pressure - center_of_mass)/diameter

def calculateTaperRatio(tip_chord, root_chord):
    return tip_chord/root_chord

def calculateAspectRatio(tip_chord, root_chord, wingspan, sweep_distance):
    wing_projected_area = (tip_chord * wingspan) + (0.5 * sweep_distance * wingspan)
    aspect_ratio = (wingspan**2) / wing_projected_area
    return aspect_ratio

def CalculateFinFlutterCriticalMachNumber(shear_modulus, aspect_ratio, fin_thickness, root_chord, coefficient, pressure, taper_ratio):
    G = shear_modulus
    AR = aspect_ratio
    P = pressure
    t = fin_thickness
    C = coefficient
    c = root_chord
    
    
    numerator = G * 2 * (AR + 2) * ((t/c)**3)
    denominator = C * (AR**3) * P * (taper_ratio + 1)
    return np.sqrt(numerator/denominator) 

def main():
    target_stability_caliber = 1.75
    center_of_pressure_from_top = calculateCenterOfPressure(target_stability_caliber, parameters.wet_COM_location_from_top, parameters.tube_outer_diameter)
    dry_stability_cal = calculateStabilityCaliber(center_of_pressure_from_top, parameters.dry_COM_location_from_top, parameters.tube_outer_diameter)

    print(f"center_of_pressure_from_top: {center_of_pressure_from_top * c.M2IN:.2f} inches from top")
    print(f"dry_stability_cal: {dry_stability_cal:.2f} cal")

    shear_modulus = 26e9 # [pa]
    fin_thickness = (1/4) * c.IN2M
    coefficient = 1.337 # what the fuck
    atmospheric_pressure = 1 * c.ATM2PA

    tip_chord = parameters.tip_chord
    root_chord = parameters.root_chord
    wingspan = parameters.wingspan
    sweep_distance = root_chord - tip_chord
    taper_ratio = calculateTaperRatio(tip_chord, root_chord)
    aspect_ratio = calculateAspectRatio(tip_chord, root_chord, wingspan, sweep_distance)

    critical_mach = CalculateFinFlutterCriticalMachNumber(shear_modulus, aspect_ratio, fin_thickness, root_chord, coefficient, atmospheric_pressure, taper_ratio)
    print(f"critical_mach: {critical_mach} Mach")


main()