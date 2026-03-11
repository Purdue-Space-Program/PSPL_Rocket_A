import numpy as np
import os
import sys
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import constants as c

#vehicle parameters
P_tank_LOx = 417 #psia #LOx tank starting pressure
P_tank_IPA = 417 #psia #IPA tank starting pressure
rho_Lox = 1141 #kg/m^3 #density of incompressible LOx
rho_IPA = 785 #kg/m^3 #density of incompressible IPA
mu_LOx = 2e-4 #kg/m/s #dynamic viscosity of incompressible LOx
mu_IPA = 0.00204 #kg/m/s #dynamic viscosity of incompressible IPA

#pipes
outer_diameter = 0.5 #inches
wall_thickness = 0.049 #inches #other available: 0.049, 0.083, 0.065 in
inner_diameter = outer_diameter - (2 * wall_thickness) #inches
inner_area = np.pi * (inner_diameter / 2) ** 2 #in^2

#mass flow rates
mdot_Lox = 1.91 * c.LBM2KG #kg/s #LOx mass flow rate
mdot_IPA = 1.91 * c.LBM2KG #kg/s #IPA mass flow rate

#line velocities 
Line_velocity_Lox = mdot_Lox / (rho_Lox * (inner_area * c.IN22M2)) #m/s
Line_velocity_IPA = mdot_IPA / (rho_IPA * (inner_area * c.IN22M2)) #m/s

#calculating Reynolds number to determine flow regime

Re_LOx = (rho_Lox * Line_velocity_Lox * (inner_diameter * c.IN2M)) / mu_LOx
Re_IPA = (rho_IPA * Line_velocity_IPA * (inner_diameter * c.IN2M)) / mu_IPA

#calculating friction factor
abs_rough = 0.0015e-3 #m absolute roughness
f_LOx = 1/((-1.8 * np.log10((((abs_rough/(inner_diameter*c.IN2M))/3.7)**1.11) + (6.9/Re_LOx)))**2) #friction factor for LOx
f_IPA = 1/((-1.8 * np.log10((((abs_rough/(inner_diameter*c.IN2M))/3.7)**1.11) + (6.9/Re_IPA)))**2) #friction factor for IPA

#length of pipes
length_LOx = ((2.5*2)+0.33) * c.IN2M #m #length of lower LOx pipe
length_IPA = 8 * c.IN2M #m #length of lower IPA pipe


#pressure drop calculation using Darcy-Weisbach equation

delta_P_LOx_pipe = (f_LOx * (length_LOx) / (inner_diameter * c.IN2M)) * (rho_Lox * Line_velocity_Lox ** 2) / 2 #Pa #pressure drop across lower LOx pipe
delta_P_IPA_pipe = (f_IPA * (length_IPA) / (inner_diameter * c.IN2M)) * (rho_IPA * Line_velocity_IPA ** 2) / 2 #Pa #pressure drop across lower IPA pipe

delta_P_LOx_loss = (0.05 + 0.4) * (rho_Lox * Line_velocity_Lox ** 2) / 2 #Pa #pressure drop across LOx losses (bends, fittings, etc.) estimated using K factors
delta_P_IPA_loss = 0 #since only a straight pipe in lower plumbing, K=0

delta_P_LOx = delta_P_LOx_pipe + delta_P_LOx_loss #Pa #total pressure drop for LOx
delta_P_IPA = delta_P_IPA_pipe + delta_P_IPA_loss #Pa #total pressure drop for IPA

print(f"LOx Line Velocity: {Line_velocity_Lox:.2f} m/s")
print(f"IPA Line Velocity: {Line_velocity_IPA:.2f} m/s")
print(f"LOx Total Pressure Drop: {delta_P_LOx:.2f} Pa")
print(f"IPA Total Pressure Drop: {delta_P_IPA:.2f} Pa")

print(f"LOx Line Velocity: {Line_velocity_Lox*c.M2FT:.2f} ft/s")
print(f"IPA Line Velocity: {Line_velocity_IPA*c.M2FT:.2f} ft/s")
print(f"LOx Total Pressure Drop: {delta_P_LOx*c.PA2PSI:.2f} psi")
print(f"IPA Total Pressure Drop: {delta_P_IPA*c.PA2PSI:.2f} psi")