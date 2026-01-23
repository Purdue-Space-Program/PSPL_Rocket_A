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

def getCoM():
    with open('vehicle_parameters.csv', mode = 'r') as file:
        reader = c.DictReader(file)
        for row in reader:
            if(row["parameter_name"] == "wet_COM_location_from_top"):
                return float(row["value"])*m_to_in

def calculateCoP():
    target_stability_caliber = 1.5
    diameter = getDiameter()
    center_of_gravity = getCoM()

    return ((target_stability_caliber)*(diameter))+center_of_gravity

def main():
    center_of_pressure = calculateCoP()
    print(center_of_pressure)

main()