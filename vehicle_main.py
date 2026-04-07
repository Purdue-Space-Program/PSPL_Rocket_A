import os
from pathlib import Path
import scipy.io as sio 
import numpy as np

import vehicle_parameters
import vehicle_parameters_functions
import print_filter
import constants as c

import SFD.six_DoF_caller as six_DoF_caller
import SFD.RDOF_2.rdof_v2 as rdof_v2
import Structures_Analysis.structural_loads as structural_loads
import shear_bolted_joints


def vehicle_analysis():
    
    repository_root_path, _ = vehicle_parameters_functions.Get_Repository_Root_Path()

    PSPL_ROCKET_A_file_path = repository_root_path / Path(f"vehicle_parameters.csv")
    six_DoF_file_path = (repository_root_path / ".." / "PSPL-6DOF"/ "TheSixDoF").resolve()
    structures_analysis_file_path = (repository_root_path / "Structures_Analysis").resolve()
    
    parameters, wet_mass_distribution, dry_mass_distribution = vehicle_parameters.main()
    parameters = six_DoF_caller.main(parameters)
    parameters = rdof_v2.main(parameters)
    parameters, wet_mass_distribution = structural_loads.main(parameters, wet_mass_distribution)
    parameters = shear_bolted_joints.main(parameters)

    vehicle_parameters_functions.ExportObjectToCSV(parameters, PSPL_ROCKET_A_file_path)
    vehicle_parameters_functions.ExportObjectToCSV(wet_mass_distribution, "wet_mass_distribution")



def main():
    vehicle_analysis()

if __name__ == "__main__":
    with print_filter.context_manager(print_everything=True, print_margins=True, print_titles=True):
        main()
