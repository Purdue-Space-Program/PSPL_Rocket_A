import numpy as np
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import constants as c

'''
Using mass flow rate: 
m_dot = mass flow rate remaining in the manifold (kg/s)
x = axial position along the manifold (m)
m_dot0 = mass flow rate entering the manifold at x=0
L = 2piR total length of the manifold (m)
V = constant velocity in the manifold (m/s)
rho = fuel density (kg/m^3)

m_dot = rho * A * V


'''

initial_area = (1 * c.IN2M) * (0.125 * c.IN2M) * np.pi #initial area of the manifold (m^2)
rho = c.DENSITY_IPA # fuel density (kg/m^3) of IPA
mdot_0 = 1.91 * c.LBM2KG # mass flow rate entering the manifold (kg/s) #initial IPA Mass flow rate

radius = 1 # manifold radius (m)
n_steps = 30 # number of angular steps around teh ring

mdot_half = mdot_0 / 2.0 # flow in each direction
velocity = mdot_half / (rho * initial_area) #constant velocity in the manifold (m/s)

dtheta = np.pi / n_steps # only half the ring

# mass flow decrease per radian (per side)
dmdot_dtheta = -mdot_half / np.pi

#array initialization
theta = np.zeros(n_steps + 1)
mdot  = np.zeros(n_steps + 1)
area  = np.zeros(n_steps + 1)

#initial conditions
theta[0] = 0.0
mdot[0]  = mdot_half
area[0]  = mdot[0] / (rho * velocity)


for i in range(n_steps):
    theta[i + 1] = theta[i] + dtheta
    mdot[i + 1]  = mdot[i] + dmdot_dtheta * dtheta

    area[i + 1] = mdot[i + 1] / (rho * velocity)

print("Step | theta (rad) | mdot (kg/s) | Area (m^2)")
print("----------------------------------------------")

for i in range(n_steps + 1):
    print(f"{i:>4} | "
          f"{theta[i]:>10.4f} | "
          f"{mdot[i]:>11.6f} | "
          f"{area[i]:>12.6e}) | "
          f"{area[i] * c.M22IN2:>12.6e}")
