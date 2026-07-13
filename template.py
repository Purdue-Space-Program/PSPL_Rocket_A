import sys
import os
from pathlib import Path

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
import Vehicle_Level.vehicle_parameters_functions as vehicle_parameters_functions
import Vehicle_Level.vehicle_parameters as vehicle_parameters
import Vehicle_Level.vehicle_main as vehicle_main
import Vehicle_Level.print_filter as print_filter


def Name_of_Script_Main_Function(parameters): # this is a separate function for reasons i dont remember but its important...
    pass

def main(parameters):
    parameters = Name_of_Script_Main_Function(parameters)
    return(parameters)

if __name__ == "__main__":
    parameters, wet_mass_distribution, dry_mass_distribution = vehicle_parameters.main()
    main(parameters)