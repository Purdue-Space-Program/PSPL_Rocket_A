import sys
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *
import vehicle_parameters as v
import csv as c

m_to_in = 39.3701

def getDiameter():
    return v.propellant_tank_outer_diameter*m_to_in

def getwetCoM():
    with open('vehicle_parameters_records/vehicle_parameters_2026-01-30_16-50-01.csv', mode = 'r') as file:
        reader = c.DictReader(file)
        for row in reader:
            if(row["parameter_name"] == "wet_COM_location_from_top"):
                return float(row["value"])*m_to_in
            
def getdryCoM():
    with open('vehicle_parameters_records/vehicle_parameters_2026-01-30_16-50-01.csv', mode = 'r') as file:
        reader = c.DictReader(file)
        for row in reader:
            if(row["parameter_name"] == "dry_COM_location_from_top"):
                return float(row["value"])*m_to_in

def calculateCoP(target_cal, CoM, diameter):
    return ((target_cal)*(diameter))+CoM

def main():
    target_stability_caliber = 2.0
    center_of_pressure_wet= calculateCoP(target_stability_caliber, getwetCoM(), getDiameter())
    center_of_pressure_dry= calculateCoP(target_stability_caliber, getdryCoM(), getDiameter())

    print(center_of_pressure_wet)
    print(center_of_pressure_dry)

main()