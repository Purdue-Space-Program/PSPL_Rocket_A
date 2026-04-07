import sys
import os
from pathlib import Path
from dataclasses import dataclass

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
import vehicle_parameters_functions
import vehicle_parameters
import vehicle_main
import print_filter


def Calculate_O_Ring_Diameter(parameters):
    
    
    @dataclass
    class ORingSeal:
        groove_width: float
        groove_height: float
        bolt_thread_size: str
        number_of_bolts: float
        shear_limit_load: float

    
    fuel_tank_o_ring_seal = ORingSeal()






def main(parameters):
    parameters = Calculate_O_Ring_Diameter(parameters)
    return(parameters)

if __name__ == "__main__":
    parameters, wet_mass_distribution, dry_mass_distribution = vehicle_parameters.main()
