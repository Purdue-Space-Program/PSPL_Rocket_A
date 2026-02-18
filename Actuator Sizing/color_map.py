import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *


def calc_actuation_time_color_plot(piston_diameter, piston_stroke_length):
    pressure = 100 * PSI2PA # Setting constant 100 psi pilot pressure
    braking_torque = 50 * LBI2NM # Setting constant 50 lb in for now, will definitely change after testing
    safety_factor = 3
    rod_mass = 3.2 * LBM2KG # Estimated from CAD, just assume constant
    shaft_diameter = 0.625 * IN2M # Just assume constant
    piston_seal_length = np.pi * piston_diameter
    shaft_seal_length = np.pi * shaft_diameter # Just assume constant
    piston_seal_area = 0.21 * IN2M * piston_seal_length # worst case scenario, 300 series
    shaft_seal_area = 0.21 * IN2M * shaft_seal_length
    
    ##### Calculating net force to see if the actuator will actuate #####
    piston_force = pressure * np.pi * ((piston_diameter**2) / 4)
    # required_torque = braking_torque * safety_factor
    arm_length = piston_stroke_length / np.sqrt(2)
    # torque = arm_length * piston_force / np.sqrt(2)
    force_valve = braking_torque * np.sqrt(2) / arm_length
    fc_piston = (4 * piston_seal_length * M2IN) * LBF2N # assuming 4, worst case, for now
    fc_shaft = (4 * shaft_seal_length * M2IN) * LBF2N
    fh_piston = (18 * piston_seal_area * M22IN2) * LBF2N # from parker oring handbook figure 5-10, lowkey using it for u-cup :skull:
    fh_shaft = (18 * shaft_seal_area * M22IN2) * LBF2N
    friction_piston = fc_piston + fh_piston
    friction_shaft = fc_shaft + fh_shaft
    f_net = piston_force - force_valve - friction_piston * 2 - friction_shaft # two seals on piston, 1 on rod
    
    ##### Actually Determining the Actuation Time #####

    if f_net <= 0 * N2LBF:
        return -1 # return -1 as actuation time as a characteristic indicator
    
    time = 0
    time_step = 0.0001
    piston_velocity = 0
    dist_travelled = 0
    valve_angle = 0

    while valve_angle <= 90:
        dist_travelled = dist_travelled + piston_velocity * time_step + 0.5 * (f_net / rod_mass) * time_step**2
        piston_velocity_new = piston_velocity + (f_net * time_step) / (rod_mass)
        piston_velocity = piston_velocity_new
        if dist_travelled == 0:
            valve_angle = 0
        else:
            valve_angle = np.degrees((np.pi / 2) - np.arctan((((arm_length / dist_travelled) - (1/np.sqrt(2))) * np.sqrt(2))))
        time += time_step
    return time 

def display_color_plot():
    bore_sizes = np.linspace(1, 3, 100) * IN2M
    strokes = np.linspace(2, 20, 100) * IN2M
    X, Y = np.meshgrid(bore_sizes, strokes)
    Z = np.zeros_like(X)
    print("This may take a while... please wait.")
    for i in range(len(strokes)):
        for j in range(len(bore_sizes)):
            Z[i, j] = calc_actuation_time_color_plot(
                bore_sizes[j],
                strokes[i]
            )
        Z[Z == -1] = np.nan  # no actuation -> gray

    cmap = plt.cm.coolwarm.copy()
    cmap.set_bad('lightgray')   # NaNs -> gray
    cmap.set_over('black')      # values above vmax -> black

    mesh = plt.pcolormesh(
        X * M2IN,
        Y * M2IN,
        Z * 1000,
        cmap=cmap,
        vmin=0,
        vmax=100  # 100 ms threshold
    )

    mesh.set_clim(0, 100)  # enforce color limits
    cbar = plt.colorbar(mesh)
    cbar.set_label("Actuation Time (ms)")
    plt.xlabel("Piston Diameter [In]")
    plt.ylabel("Stroke Length [In]")
    plt.title("Actuation Time")
    plt.tight_layout()
    plt.show()

display_color_plot()