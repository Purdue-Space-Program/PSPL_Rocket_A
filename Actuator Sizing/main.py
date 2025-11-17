import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *

def calc_torque_piston(braking_torque, safety_factor, piston_force_at_200psi, piston_stroke_length):
    required_torque = braking_torque * safety_factor
    armlength = piston_stroke_length / np.sqrt(2)
    torque = armlength * piston_force_at_200psi / np.sqrt(2)
    print(f"The piston will produce ~{torque:.2f} in-lb torque at 200 psi.")
    print(f"The required torque with a safety factor of 3 is {required_torque} in-lb.")
    print(f"Length of valve arm would be {armlength} in.")
    return required_torque, armlength, torque

def actuation_time(armlength, braking_torque, torque, piston_mass):
    net_torque = torque - braking_torque
    F_net = net_torque * np.sqrt(2) / armlength
    time = 0
    time_step = 0.0001
    piston_velocity = 0
    dist_travelled = 0
    valve_angle = 0
    time_history = []
    angle_history = []
    velocity_history = []
    """while valve_angle <= 90:
        piston_velocity_new = piston_velocity + (F_net * time_step) / (piston_mass)
        dist_travelled += ((piston_velocity_new + piston_velocity) / 2) * time_step
        piston_velocity = piston_velocity_new
        time += time_step
        valve_angle = np.degrees(np.arctan(1 / (((armlength / dist_travelled) + (1/np.sqrt(2))) * np.sqrt(2))))
        time_history.append(time)
        angle_history.append(valve_angle)"""
    
    while valve_angle <= 90:
        dist_travelled = piston_velocity * time + 0.5 * (F_net / piston_mass) * time**2
        piston_velocity_new = piston_velocity + (F_net * time_step) / (piston_mass)
        piston_velocity = piston_velocity_new
        valve_angle = np.degrees((np.pi / 2) - np.arctan((((armlength / dist_travelled) - (1/np.sqrt(2))) * np.sqrt(2))))
        time += time_step
        time_history.append(time)
        angle_history.append(valve_angle)
        velocity_history.append(piston_velocity)

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
    plt.show()

    print(f"Actuation time: {time}")
    return None
# Shortlisted Piston: https://www.mcmaster.com/6498K297/
braking_torque = 240
safety_factor = 3
piston_force_at_200psi = 620
piston_stroke_length = 2.5
required_torque, armlength, torque = calc_torque_piston(braking_torque, safety_factor, piston_force_at_200psi, piston_stroke_length)
piston_mass = 5 # # [lb]
actuation_time(armlength, braking_torque, torque, piston_mass)