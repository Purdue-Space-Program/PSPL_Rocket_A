import sys
import os
from pathlib import Path

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
import vehicle_parameters_functions
import vehicle_parameters
import vehicle_main
import print_filter


def name_of_script_main_function(parameters):
    pass

def main(parameters):
    parameters = name_of_script_main_function(parameters)
    return(parameters)

if __name__ == "__main__":
    parameters, wet_mass_distribution, dry_mass_distribution = vehicle_parameters.main()
