### This code was adapted from the sim made by Hugo Filmer to work for Rocket A. David Gustafsson is the one who adapted this

# This code simulates the pressurization process in the propellant tanks over time, assuming a constant pressure in each tank.
# Its primary use is to determine the final COPV pressure and therefore if the COPV is large enough.
# Hugo Filmer

from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
import traceback

from heat_transfer_functions import *

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
import vehicle_parameters

def main(property_to_use = "density"):
    parameters, wet_mass_distribution, dry_mass_distribution = vehicle_parameters.main()

    # Simulation settings
    T_AMBIENT = 293 # [K] ambient temperature
    LOITER_TIME = 0 # [s] time between pre-pressurization and the start of flow
    LAG_TIME = 0 # [s] time the simulation should continue to run for after the run valves are closed
    DT = 0.001 # [s] simulation step size
    Q_FACTOR = 2 # [] factor to multiply heat transfer by (for conservatism)
    TEXT_OUTPUT = True # True to print summary output, including conservation and EoS checks
    PLOT_OUTPUT = True # True to make pretty plots of the results :)

