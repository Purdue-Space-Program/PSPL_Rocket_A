"""
Output Formatter
Author: Aroldo Lugo Quintanilla

Description: Takes arrays and other data and display the outputs at the sensor points
"""
import numpy as np
import pandas as pd

# Conversion Constants
LBF2N = 1 / 0.224809
FT2M = 1 / 3.28084
IN2M = 0.0254
LB2KG = 0.453592

# Terminal Text Formatting Constants
WARNING = '\033[93m'
BOLD = '\033[1m'
ENDC = '\033[0m'
OKGREEN = '\033[92m'

def outputFormatter(shear, bending, axial, shock_load, mass, length, rocket_dict_aroldo, sheet_name="Recovery Forces"):

    # Convert Units to Imperial for Outputs
    mass /= LB2KG # kg to lbm
    length_in = length/IN2M # m to in
    shear_array = shear / LBF2N # N to lbf
    bending_array =  bending / (LBF2N * FT2M) # N*m to ft-lbs
    axial_array = axial / LBF2N # N to lbf

    max_bending = np.max(bending_array)
    max_shear = np.max(shear_array)

    slice_length = (5/1000)/IN2M
    array_length = len(shear_array)
    array2 = [slice_length] * array_length
    lengths = np.cumsum(array2)

    print("-------------------------------")
    print(f"{BOLD}{sheet_name.upper} Parachute Outputs:{ENDC}\n")

    print(f"{BOLD}Total Mass:{ENDC} {mass:.2f} lbm")
    print(f"{BOLD}Total Length:{ENDC} {length_in:.2f} in\n")
    print(f"{BOLD}Max Shear Force:{ENDC} {max_shear:.2f} lbf")
    print(f"{BOLD}Max Bending Moment:{ENDC} {max_bending:.2f} ft-lbs")
    if shock_load: print(f"{BOLD}Main Axial Shock Load:{ENDC} {shock_load/LBF2N:.2f} lbf\n")
    # print(f"{BOLD}Max Shear Shock Load:{ENDC} {shock_load[1]:.2f} lbf\n")

    sensor_labels = ["Recovery Bay-Helium Bay", "Helium Bay-Upper Airframe", "Upper Airframe-Lox Tank", 
        "Lox Tank-Mid Airframe", "Mid Airframe-Fuel Tank", "Fuel Tank-Lower Airframe", "Lower Airframe-Engine", "Aft of Rocket"
    ] #"Recovery Bay-Helium Bay",

    totalLength = length
    index = []
    x = 0

    sensor_df = pd.DataFrame({})

    for i, (key, section) in enumerate(rocket_dict_aroldo.items()):
        # if sensor_labels[x] == 'Aft of Rocket':
        #     break
        print(key, totalLength, section['length'])
        if x >= len(sensor_labels): break
        if 'length'==0: continue
        if key == 'recovery_bay': continue
        totalLength -= section['length']
        forces_index = int(totalLength / .005)
        if forces_index < 0:
            forces_index = 0
        axial_force = axial_array[forces_index]
        shear_force = shear_array[forces_index]
        bending_moment = bending_array[forces_index] 
        # bottom_peq, max_peq = peq[i]
        label = sensor_labels[x]

        index.append(forces_index) #Finding the slice length
        print("-------------------------------")
        print(f'{BOLD}{sensor_labels[x]}{ENDC}')
        print(f"Distance from aft: {lengths[forces_index]:.2f} in")
        print(f"Axial Compression: {axial_force:.2f} lbf")
        print(f"Bending Moment: {bending_moment:.2f} ft-lbs")
        print(f"Shear Force: {shear_force:.2f} lbf")
        sensor_df.loc[label, 'Axial (lbf)'] = axial_force
        sensor_df.loc[label, 'Shear (lbf)'] = shear_force
        sensor_df.loc[label, 'Bending Moment (ft-lbs)'] = bending_moment
        # sensor_df.loc[label, 'Bottom PEQ'] = bottom_peq
        # sensor_df.loc[label, 'Max PEQ'] = bottom_peq
        x += 1


    # Save to An Excel File --> outputs.xslx
    data = pd.DataFrame({"Distance from Aft (in)": np.round(lengths, decimals=2), "Shear (lbf)": np.round(shear_array, decimals=2), "Bending (ft-lbs)": np.round(bending_array, decimals=2), "Axial (lbf)": np.round(axial_array, decimals=2)})
    print()

    try:
        with pd.ExcelWriter("outputs.xlsx", mode='a', if_sheet_exists='replace') as outputs:
            data.to_excel(outputs, sheet_name=sheet_name, index=False)
            sensor_df.to_excel(outputs, sheet_name=sheet_name + " - Sensor Points")

        print(f'{OKGREEN}{sheet_name} Outputs saved to "outputs.xlsx"{ENDC}')
    except PermissionError:
        print(f'{WARNING}"outputs.xlsx" is being used by another application. Please close it, then run RDOF again.{ENDC}')
