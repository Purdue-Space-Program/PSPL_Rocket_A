import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *

###############################
# CHOOSE MODE "test" or "real"
piston = "test" # test or real
###############################

###############################
# INPUT PARAMETERS
###############################
braking_torque = 240 * LBI2NM 
safety_factor = 3
piston_stroke_length = 2 * IN2M 
piston_mass = 1.6 * LBM2KG # Estimated from CAD
piston_diameter = 17/16 * IN2M
piston_retracted_length = 9.44 * IN2M
piston_extended_length = piston_retracted_length + piston_stroke_length
shaft_diameter = 0.312 * IN2M
piston_seal_length = np.pi * piston_diameter
shaft_seal_length = np.pi * shaft_diameter
piston_seal_area = 0.21 * IN2M * piston_seal_length # worst case scenario, 300 series 
shaft_seal_area = 0.21 * IN2M * shaft_seal_length
pressure = 200 * PSI2PA
piston_force = pressure * np.pi * piston_diameter**2 / 4
print(f'Maximum possible net force disregarding friction (and valve arm if real condition): {piston_force * N2LBF}')

def calc_net_force_real(piston_force, piston_seal_length, shaft_seal_length, piston_seal_area, shaft_seal_area, braking_torque, armlength):
    fc_piston = (4 * piston_seal_length * M2IN) * LBF2N # assuming 4, worst case, for now
    fc_shaft = (4 * shaft_seal_length * M2IN) * LBF2N
    fh_piston = (18 * piston_seal_area * M22IN2) * LBF2N # from parker oring handbook figure 5-10, lowkey using it for u-cup :skull:
    fh_shaft = (18 * shaft_seal_area * M22IN2) * LBF2N
    friction_piston = fc_piston + fh_piston
    friction_shaft = fc_shaft + fh_shaft
    force_valve = braking_torque * np.sqrt(2) / armlength
    print(f"fc_piston: {fc_piston * N2LBF} LBF")
    print(f"fc_shaft: {fc_shaft * N2LBF} LBF")
    print(f"fh_piston: {fh_piston * N2LBF} LBF")
    print(f"fh_shaft: {fh_shaft * N2LBF} LBF")
    f_net = piston_force - force_valve - friction_piston * 2 - friction_shaft # two seals on piston, 1 on rod
    print(f"F_net: {f_net * N2LBF} LBF")
    return f_net

def calc_net_force_test(piston_force, piston_seal_length, shaft_seal_length, piston_seal_area, shaft_seal_area):
    fc_piston = (4 * piston_seal_length * M2IN) * LBF2N # assuming 4, worst case, for now
    fc_shaft = (4 * shaft_seal_length * M2IN) * LBF2N
    fh_piston = (18 * piston_seal_area * M22IN2) * LBF2N # from parker oring handbook figure 5-10, lowkey using it for u-cup :skull:
    fh_shaft = (18 * shaft_seal_area * M22IN2) * LBF2N
    friction_piston = fc_piston + fh_piston
    friction_shaft = fc_shaft + fh_shaft
    print(f"fc_piston: {fc_piston * N2LBF} LBF")
    print(f"fc_shaft: {fc_shaft * N2LBF} LBF")
    print(f"fh_piston: {fh_piston * N2LBF} LBF")
    print(f"fh_shaft: {fh_shaft * N2LBF} LBF")
    f_net = piston_force - friction_piston * 2 - friction_shaft # two seals on piston, 1 on rod
    print(f"F_net: {f_net * N2LBF} LBF")
    return f_net

def calc_volumetric_flow(volume_swept_history, time_history):
    volumetric_flow_history = []
    time_step = time_history[1] - time_history[0] if len(time_history) > 1 else 0
    volumetric_flow_history.append(0)
    for i in range(1, len(volume_swept_history)):
        volumetric_flow = (volume_swept_history[i] - volume_swept_history[i-1]) / time_step
        volumetric_flow_history.append(volumetric_flow)
    plt.subplot(2, 1, 1)
    plt.plot(volume_swept_history, volumetric_flow_history)
    plt.title("Volume vs Volumetric Flow")
    plt.subplot(2, 1, 2)
    plt.plot(time_history, volumetric_flow_history)
    plt.title("Time vs Volumetric Flow")
    plt.tight_layout()
    plt.show()
    return volumetric_flow_history, time_history

def calc_torque_piston(braking_torque, safety_factor, piston_force, piston_stroke_length):
    required_torque = braking_torque * safety_factor
    armlength = piston_stroke_length / np.sqrt(2)
    torque = armlength * piston_force / np.sqrt(2)
    print(f"The piston will produce ~{torque * NM2LBI} torque at 200 psi.")
    print(f"The required torque with a safety factor of 3 is {required_torque * NM2LBI}")
    print(f"Length of valve arm would be {armlength * M2IN}")
    return required_torque, armlength, torque

def actuation_time_vol_flow(piston_diameter, volumetric_flow_history, time_history, piston_stroke_length):
    time = 0
    time_history_vol_flow = []
    distance_travelled_history = []
    distance_travelled = 0
    piston_area = np.pi * piston_diameter**2 / 4 
    time_step = time_history[1] - time_history[0]
    for volumetric_flow in volumetric_flow_history:
        piston_velocity = volumetric_flow / piston_area if volumetric_flow is not None else 0
        distance_travelled = distance_travelled + piston_velocity * time_step
        distance_travelled_history.append(distance_travelled)
        time_history_vol_flow.append(time)
        time = time + time_step
        if distance_travelled >= piston_stroke_length:
            return distance_travelled, time 
    else:
        print('Ran out of volumetric flow values in volumetric flow history, something is broken :(')
        return 0, 0
    
def actuation_time_kinematics_real(F_net, piston_mass, piston_diameter, armlength):
    time = 0 
    time_step = 0.00001  
    piston_velocity = 0 
    dist_travelled = 0
    valve_angle = 0
    volume_swept = 0
    time_history = []
    angle_history = []
    velocity_history = []
    volume_swept_history = []
    distance_travelled_history = []
    while valve_angle <= 90:
        dist_travelled = dist_travelled + piston_velocity * time_step + 0.5 * (F_net / piston_mass) * time_step**2
        piston_velocity_new = piston_velocity + (F_net * time_step) / (piston_mass)
        piston_velocity = piston_velocity_new
        if dist_travelled == 0:
            valve_angle = 0
        else:
            valve_angle = np.degrees((np.pi / 2) - np.arctan((((armlength / dist_travelled) - (1/np.sqrt(2))) * np.sqrt(2))))
        volume_swept = dist_travelled * np.pi * (piston_diameter / 2)**2
        time_history.append(time)
        angle_history.append(valve_angle)
        velocity_history.append(piston_velocity)
        volume_swept_history.append(volume_swept)
        distance_travelled_history.append(dist_travelled)
        time += time_step
    plt.subplot(1, 2, 1)
    plt.plot(time_history, angle_history)
    plt.xlabel("Actuation Time")
    plt.ylabel("Valve Angle")
    plt.title("Valve Angle Over Time")
    plt.ylim(0, 90)
    plt.xlim(0, time)

    plt.subplot(2, 2, 2)
    plt.plot(time_history, velocity_history)
    plt.xlabel("Time")
    plt.ylabel("Velocity")
    plt.title("Time vs Velocity")

    plt.subplot(2, 2, 3)
    plt.plot(time_history, volume_swept_history)
    plt.xlabel("Time")
    plt.ylabel("Volume")
    plt.title("Time vs Volume Swept")

    plt.subplot(2, 2, 4)
    plt.plot(time_history, distance_travelled_history)
    plt.xlabel("Time")
    plt.ylabel("Distance Swept")
    plt.title("Time vs Distance Swept")
    plt.tight_layout()
    plt.show()

    print(f"Actuation time: {time}s")
    print(f"Stroke length when using valve angle condition: {distance_travelled_history[-1] * M2IN} in")

    return volume_swept_history, time_history

def actuation_time_kinematics_test(F_net, piston_mass, piston_diameter, piston_stroke_length):
    time = 0 
    time_step = 0.00001  
    piston_velocity = 0 
    dist_travelled = 0
    volume_swept = 0
    time_history = []
    velocity_history = []
    volume_swept_history = []
    distance_travelled_history = []
    while dist_travelled <= piston_stroke_length:
        dist_travelled = dist_travelled + piston_velocity * time_step + 0.5 * (F_net / piston_mass) * time_step**2
        piston_velocity_new = piston_velocity + (F_net * time_step) / (piston_mass)
        piston_velocity = piston_velocity_new
        volume_swept = dist_travelled * np.pi * (piston_diameter / 2)**2
        time_history.append(time)
        velocity_history.append(piston_velocity)
        volume_swept_history.append(volume_swept)
        distance_travelled_history.append(dist_travelled)
        time += time_step

    plt.subplot(1, 3, 1)
    plt.plot(time_history, velocity_history)
    plt.xlabel("Time")
    plt.ylabel("Velocity")
    plt.title("Time vs Velocity")

    plt.subplot(1, 3, 2)
    plt.plot(time_history, volume_swept_history)
    plt.xlabel("Time")
    plt.ylabel("Volume")
    plt.title("Time vs Volume Swept")

    plt.subplot(1, 3, 3)
    plt.plot(time_history, distance_travelled_history)
    plt.xlabel("Time")
    plt.ylabel("Distance Swept")
    plt.title("Time vs Distance Swept")
    plt.tight_layout()
    plt.show()

    print(f"Actuation time: {time}s")
    print(f"Stroke length when using valve angle condition: {distance_travelled_history[-1] * M2IN} in")

    return volume_swept_history, time_history

# Shortlisted Piston: https://pspliquids.slack.com/archives/C09C5J1EJDB/p1764894397354269?thread_ts=1764888234.600949&cid=C09C5J1EJDB

if piston.lower() == "test":
    f_net = calc_net_force_test(piston_force, piston_seal_length, shaft_seal_length, piston_seal_area, shaft_seal_area)
    volume_swept_history, time_history = actuation_time_kinematics_test(f_net, piston_mass, piston_diameter, piston_stroke_length)
elif piston.lower() == "real":
    required_torque, armlength, torque = calc_torque_piston(braking_torque, safety_factor, piston_force, piston_stroke_length)
    f_net = calc_net_force_real(piston_force, piston_seal_length, shaft_seal_length, piston_seal_area, shaft_seal_area, braking_torque, armlength)
    volume_swept_history, time_history = actuation_time_kinematics_real(f_net, piston_mass, piston_diameter, armlength)
else:
    print('Invalid piston chosen')
volumetric_flow_history, time_history = calc_volumetric_flow(volume_swept_history, time_history)
distance_travelled, time = actuation_time_vol_flow(piston_diameter, volumetric_flow_history, time_history, piston_stroke_length)
print(f"Actuation time using vol flow: {time}s")