import numpy as np
import math as m
import matplotlib.pyplot as plt
from vehicle_pressurization_simulation import crit_length_inputs
import CoolProp.CoolProp as cool


gamma = 1.47
R = 8314 / 28.014 # Specific gas constant for N2 [J/kgK]
e = 1.5e-6 # Surface roughness [m] check again
# Inlet Properties         
P_static_pa, T_static_k, mdot_lox, mdot_IPA, times = crit_length_inputs()
mdot_before_reg = mdot_lox + mdot_IPA

def dyn_viscosity(T,P):
    # T in Kelvin, returns viscosity in Pa-s
    dyn_viscosity = cool.PropsSI('V', 'T', T, 'P', P, 'Nitrogen')
    return dyn_viscosity

def flow_velocity(m_dot, inner_dia, density):
    """Calculates flow velocity based on mass flow and density."""
    area = (m.pi * inner_dia**2) / 4
    return (m_dot / (density * area))

def mach_at_point(velocity, gamma, R, T_static):
    """Calculates Mach number: M = V / sqrt(gamma * R * T)."""
    speed_of_sound = np.sqrt(gamma * R * T_static)
    return velocity / speed_of_sound

def f_lstar_ratio(mach, gamma):
    f_lstar_ratio = ((1-mach**2)/(gamma*mach**2)) + ((gamma+1)/(2*gamma))*np.log(((gamma+1)*mach**2)/ (2 + (gamma-1)*mach**2))
    return f_lstar_ratio

def reynolds_num(density, velocity, dynamic_viscosity, inner_dia):
    """Calculates Reynolds number."""
    return (density * velocity * inner_dia) / dynamic_viscosity

def fric_factor(re_num, e, inner_dia):
    """Haaland equation for Darcy Friction Factor."""
    # Handle potential division by zero or very low Re
    term = ((e / inner_dia) / 3.7)**1.11 + (6.9 / re_num)
    return (-1.8 * np.log10(term))**-2

# Geometry: Inner diameters converted to meters
# inner_diameters_m = 10**(-2) * np.array([
#     0.4572, 0.7036, 1.0211, 0.9398, 0.8484,
#     1.2573, 1.5748, 1.4224, 1.2954, 2.2098,
#     2.1184, 2.0574, 1.9304, 1.5850, 3.2004])
inner_diameters_m = 10**(-2) * np.array([
    1.0211, 0.9398, 0.8484])

#in_dia = np.linspace(0.4, 3.2, 10)

for dia in inner_diameters_m:
    inner_dia = dia
    R = 8314 / 28.014 
    gamma = 1.47 
    P_static, T_static_kel, mdot_lox, mdot_IPA, times =  crit_length_inputs()
    mdot_before_reg = mdot_lox + mdot_IPA
    # T_static_k = T_static_kel
    # P_static_pa = P_static
    rho_nitrogen = P_static_pa / (R * T_static_k) # Ideal Gas Law
    visc = cool.PropsSI('viscosity', 'T', 300, 'P', 2.8958e7, 'Nitrogen')
    vel_total =  flow_velocity(mdot_before_reg, inner_dia, rho_nitrogen)
    mach_total = mach_at_point(vel_total, gamma, R, T_static_k)
    re_lox = reynolds_num(rho_nitrogen, vel_total, visc, inner_dia)
    f_beforereg = fric_factor(re_lox, e, inner_dia)
    f_lstar_ratio_beforereg = f_lstar_ratio(mach_total, gamma)
    L_star_before_reg = (f_lstar_ratio_beforereg * inner_dia) / f_beforereg
    dia_label = f'Dia: {dia *39.37:.3f} in'
    plt.plot(times,L_star_before_reg*1e-3*3.28084, marker='o')
    #plt.legend(loc='best')


plt.title('Critical Pipe Length (COPV to Regulator) vs. Time', fontsize=14)
plt.xlabel('Time (s)', fontweight='bold')
plt.ylabel('Critical Pipe Length (COPV to Regulator) (ft)', fontweight='bold')
plt.grid(True, which="both", linestyle='--', alpha=0.5)
plt.show()


enthalpy = cool.PropsSI('H', 'T', T_static_k, 'P', P_static_pa, 'Nitrogen')
print(enthalpy)






# fig1, ax1 = plt.subplots(figsize=(10, 6))
# #fig2, ax2 = plt.subplots(figsize=(10, 6))
# for dia in inner_diameters_m: 
#     inner_dia = dia
#     # Inlet Properties  
#     R = 8314 / 28.014
#     gamma = 1.47 
#     P_static, T_static_kel, mdot_lox, mdot_IPA, times =  crit_length_inputs()
#     T_static_k = T_static_kel
#     P_static_pa = P_static
#     rho_nitrogen = P_static_pa / (R * T_static_k)# Ideal Gas Law
#     # Perform Calculations
#     # visc = np.zeros(240)
#     # T_static = np.array([T_static_k])
#     # P_static = np.array([P_static_pa])
#     # print(T_static)
#     # for i in len(240):
#     #     visc(i-1) == dyn_viscosity(T_static(i), P_static(i))
#     print(T_static_k)
#     print(P_static_pa)
#     visc = 1.78e-5 
        
#     #print (rho_nitrogen)
#     # Velocity and Mach
#      ## calculate flow velocities for every inner diameter, create a for loop for all the things to be calculated new time in loop
#     vel_lox = flow_velocity(mdot_lox, inner_dia, rho_nitrogen)
#     vel_fuel = flow_velocity(mdot_IPA, inner_dia, rho_nitrogen)
#     # mach_lox = mach_at_point(vel_lox, gamma, R, T_static_k)
#     # print(mach_lox)

#      # Reynolds and Friction
#     re_lox = reynolds_num(rho_nitrogen, vel_lox, visc, inner_dia)
#     f_lox = fric_factor(re_lox, e, inner_dia)
#     #length_lox = 2*P_static_pa*inner_dia / (rho_nitrogen*f_lox*vel_lox**2)

#     re_fuel = reynolds_num(rho_nitrogen, vel_fuel, visc, inner_dia)
#     f_fuel = fric_factor(re_fuel, e, inner_dia)
#     #length_fuel = 2*P_static_pa*inner_dia / (rho_nitrogen*f_fuel*vel_fuel**2)
#      # Calculate L_star (Maximum allowable pipe length before choking)
#     # L_star = (Ratio * D) / f
#     f_lstar_ratio_lox = f_lstar_ratio(mach_lox, gamma)
#     L_star_lox = (f_lstar_ratio_lox * inner_dia) / f_lox
#     #print(L_star_lox)   
#     plt.plot(times,length_lox, '.')
#     plt.legend(loc='best')

#     # plt.plot(times,length_fuel*3.28084, '.')
#     # # plt.legend(loc='best')
   
    

# # 3. Move formatting and show OUTSIDE the loop
# # plt.title('Critical Pipe Length (COPV to Lox tank) vs. Time', fontsize=14)
# # plt.xlabel('Time (s)', fontweight='bold')
# # plt.ylabel('Critical Pipe Length (COPV to Lox tank) (m)', fontweight='bold')
# # #plt.yscale('log') 
# # plt.grid(True, which="both", linestyle='--', alpha=0.5)
# # plt.show()

# plt.title('Critical Pipe Length (COPV to IPA tank) vs. Time', fontsize=14)
# plt.xlabel('Time (s)', fontweight='bold')
# plt.ylabel('Critical Pipe Length (COPV to IPA tank) (ft)', fontweight='bold')
# #plt.yscale('log') 
# plt.grid(True, which="both", linestyle='--', alpha=0.5)
# plt.show()

