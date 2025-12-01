import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *

def calc_net_force():
    return None

def calc_volumetric_flow(volume_swept_history, time_history):
    volumetric_flow_history = []
    for volume, time in zip(volume_swept_history, time_history):
        volumetric_flow = volume / time
        volumetric_flow_history.append(volumetric_flow)
    plt.subplot(2, 1, 1)
    plt.plot(volume_swept_history, volumetric_flow_history)
    plt.title("Volume vs Volumetric Flow")
    plt.subplot(2, 1, 2)
    plt.plot(time_history, volumetric_flow_history)
    plt.title("Time vs Volumetric Flow")
    plt.tight_layout()
    plt.show()
    return None

def calc_torque_piston(braking_torque, safety_factor, piston_force_at_200psi, piston_stroke_length):
    required_torque = braking_torque * safety_factor
    armlength = piston_stroke_length / np.sqrt(2)
    torque = armlength * piston_force_at_200psi / np.sqrt(2)
    print(f"The piston will produce ~{torque * NM2LBI} torque at 200 psi.")
    print(f"The required torque with a safety factor of 3 is {required_torque * NM2LBI}")
    print(f"Length of valve arm would be {armlength * M2IN}")
    return required_torque, armlength, torque

def actuation_time(armlength, braking_torque, torque, piston_mass, piston_diameter):
    net_torque = torque - braking_torque
    F_net = net_torque * np.sqrt(2) / armlength
    time = 0 
    time_step = 0.0001 
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
        dist_travelled = piston_velocity * time + 0.5 * (F_net / piston_mass) * time**2
        piston_velocity_new = piston_velocity + (F_net * time_step) / (piston_mass)
        piston_velocity = piston_velocity_new
        if dist_travelled ==0:
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

    plt.subplot(2, 2, 1)
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

    print(f"Actuation time: {time}")
    return volume_swept_history, time_history

# Shortlisted Piston: https://www.mcmaster.com/6498K297/

breaking_torque = 240 * LBI2NM 
safety_factor = 3
piston_force_at_200psi = 620 * LBF2N 
piston_stroke_length = 2.5 * IN2M 
required_torque, armlength, torque = calc_torque_piston(breaking_torque, safety_factor, piston_force_at_200psi, piston_stroke_length)
piston_mass = 5 * LBM2KG 
piston_diameter = 2 * IN2M
volume_swept_history, time_history = actuation_time(armlength, breaking_torque, torque, piston_mass, piston_diameter)
calc_volumetric_flow(volume_swept_history, time_history)