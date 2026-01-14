import numpy as np

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

mdot_0 = 1 # mass flow rate entering the manifold (kg/s)
rho = 1 # fuel density (kg/m^3)
velocity = 5.0 # constant velocity (m/s)
radius = 1 # manifold radius (m)
n_steps = 20 # number of angular steps around teh ring

circumference = 2.0 * np.pi * radius
dtheta = 2.0 * np.pi / n_steps

'''
mass flow rate decreases linearly with angle theta:
Over the full circle, the mass flow rate decreases from mdot_0 to 0:
'''

# mass flow loss per radian
dmdot_dtheta = -mdot_0 / (2.0 * np.pi)

# arrays to store values at each step
theta = np.zeros(n_steps + 1)
mdot = np.zeros(n_steps + 1)
area = np.zeros(n_steps + 1)

# initial conditions
theta[0] = 0
mdot[0] = mdot_0
area[0] = mdot[0] / (rho * velocity)

'''
Euler's method:

increase angle
reduce mass flow rate by dmdot_dtheta
new area calculate
'''
for i in range(n_steps):
    theta[i+1] = theta[i] + dtheta
    mdot[i+1]  = mdot[i] + (dmdot_dtheta * dtheta)
    area[i+1]  = mdot[i+1] / (rho * velocity)

# output
print("Step | theta (rad) | mdot (kg/s) | Area (m^2)")
print("----------------------------------------------")

for i in range(n_steps + 1):
    print(f"{i:>4} | "
          f"{theta[i]:>10.4f} | "
          f"{mdot[i]:>11.4f} | "
          f"{area[i]:>10.6e}")
