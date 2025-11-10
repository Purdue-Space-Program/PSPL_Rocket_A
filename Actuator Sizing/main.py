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
    return None

braking_torque = 240
safety_factor = 3
piston_force_at_200psi = 620
piston_stroke_length = 2.5
calc_torque_piston(braking_torque, safety_factor, piston_force_at_200psi, piston_stroke_length)