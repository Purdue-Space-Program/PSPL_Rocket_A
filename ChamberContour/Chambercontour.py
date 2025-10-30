import numpy as np
import matplotlib.pyplot as plt

def nozzle_contour(Dt, exp_ratio, Lstar, contract_ratio, con_angle, Dc, filename):
    # Unit conversions
    IN2M = 0.0254
    
    # Nozzle dimensions
    Lstar = Lstar * IN2M
    Dc = Dc * IN2M
    Dt = Dt * IN2M
    Rt = Dt/2  # Throat radius
    At = np.pi*Rt**2
    Ae = At * exp_ratio # Exit area
    Re = np.sqrt(Ae / np.pi)  # Exit radius  
    Ac = np.pi/4 * Dc**2
    Rc = Dc / 2
    percent_rao = 0.8
    Ln = percent_rao * ((np.sqrt(exp_ratio) - 1) * Rt) / np.tan(np.deg2rad(15)) # Nozzle length
    theta_n = np.deg2rad(-0.00001*exp_ratio**4 + 0.0017*exp_ratio**3 - 0.0775*exp_ratio**2 + 1.5769*exp_ratio + 16.805) # Initial bell angle
    theta_e = np.deg2rad(0.00001*exp_ratio**4 - 0.0015*exp_ratio**3 + 0.0665*exp_ratio**2 - 1.2849*exp_ratio + 17.895) # Final bell angle

    if theta_e < 0:
        print("Exit angle is negative")
        theta_e = 0
# http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
# Throat section
    # Throat converging section arc
    con_angle = np.deg2rad(con_angle)
    theta_con = np.linspace(-np.pi/2 - con_angle, -np.pi/2)
    x_con = 1.5 * Rt * np.cos(theta_con) # circle with radius 1.5Rt//
    y_con = (1.5 * Rt * np.sin(theta_con)) + (1.5 * Rt) + Rt
    x_con = np.delete(x_con, -1)
    y_con = np.delete(y_con, -1)
    
    # Throat diverging section arc
    theta_div = np.linspace(-np.pi/2, theta_n - np.pi/2)
    x_div = 0.382 * Rt * np.cos(theta_div) #circle with radius 0.382Rt//
    y_div = (0.382 * Rt * np.sin(theta_div)) + (0.382 * Rt) + Rt
    x_div = np.delete(x_div, -1)
    y_div = np.delete(y_div, -1)
    
# Bell section(diverging)
    Nx = 0.382 * Rt * np.cos(theta_n - np.pi/2)
    Ny = (0.382 * Rt * np.sin(theta_n - np.pi/2)) + (0.382 * Rt) + Rt
    Ex = Ln
    Ey = Re
    m1 = np.tan(theta_n)
    m2 = np.tan(theta_e)
    c1 = Ny - (m1 * Nx)
    c2 = Ey - (m2 * Ex)
    Qx = (c2 - c1) / (m1 - m2)
    Qy = ((m1 * c2) - (m2 * c1)) / (m1 - m2)
    
    theta_bell = np.linspace(0, 1)
    x_bell = (((1 - theta_bell)**2) * Nx) + (2 * (1 - theta_bell) * theta_bell * Qx) + ((theta_bell**2) * Ex)
    y_bell = (((1 - theta_bell)**2) * Ny) + (2 * (1 - theta_bell) * theta_bell * Qy) + ((theta_bell**2) * Ey) #bezier curve to draw parabola
# Converging section
    # Converging cone
    Rconv = 1 * Rc
    y_length_cone = Rc - Rconv * (np.sin(np.pi/2) - np.sin(np.pi/2 - con_angle)) - y_con[0] #//
    x_length_cone = (y_length_cone) / np.tan(-con_angle)
    x_cone = np.linspace(x_length_cone, 0) + x_con[0]
    y_cone = np.linspace(y_length_cone, 0) + y_con[0]
    x_cone = np.delete(x_cone, -1)
    y_cone = np.delete(y_cone, -1)
    
    # Converging radius
    theta_con_rad = np.linspace(np.pi/2, np.pi/2 - con_angle)
    x_con_rad = Rconv * np.cos(theta_con_rad) - Rconv * np.cos(np.pi/2 - con_angle) + x_cone[0]
    y_con_rad = Rconv * np.sin(theta_con_rad) - Rconv + Rc
    x_con_rad = np.delete(x_con_rad, -1)
    y_con_rad = np.delete(y_con_rad, -1)

# Chamber section
    # Volume calculations
    V_con = np.pi * np.trapezoid(np.concatenate([y_con_rad, y_cone, y_con])**2, x=np.concatenate([x_con_rad, x_cone, x_con])) # Volume of converging section
    Vc = (Lstar * At) - V_con # Volume of cylindrical section 
    Lc = Vc / Ac # Length of chamber
    # Chamber
    x_c = np.linspace(-Lc+x_con_rad[0], x_con_rad[0], num=3) 
    y_c = np.zeros(np.size(x_c)) + Rc
    x_c = np.delete(x_c, -1)
    y_c = np.delete(y_c, -1)

    x_arr = np.concatenate([x_c, x_con_rad, x_cone, x_con, x_div, x_bell])
    y_arr = np.concatenate([y_c, y_con_rad, y_cone, y_con, y_div, y_bell])

    first_point_depth = x_arr[0]
    last_point_depth = x_arr[-1]
    chamber_length = last_point_depth - first_point_depth

    print(f'Combustion chamber length: {Lc / IN2M:.3f} in')
    print(f'Total chamber length: {chamber_length / IN2M:.3f} in')
    print(f'Injector-to-throat length: {x_arr[0]*-1 / IN2M:.3f} in')
    print(f'Initial parabola angle: {np.rad2deg(theta_n):.2f} degrees')
    print(f'Exit parabola angle: {np.rad2deg(theta_e):.2f} degrees')
    
    z_arr = np.zeros(np.size(x_arr))
    nozzle = np.transpose(np.array([x_arr, y_arr, z_arr]))
    nozzle_in = nozzle / IN2M
    np.savetxt(filename, nozzle_in, delimiter=',')


# Nozzle contour plot
    plt.plot(x_arr/IN2M, y_arr/IN2M, 'k', x_arr/IN2M, -y_arr/IN2M, 'k')
    
    print(f"Start of exit parabola: ({Nx/ IN2M:.4f}, {Ny/ IN2M:.4f})")
    print(f"Middle of exit parabola: ({Qx/ IN2M:.4f}, {Qy/ IN2M:.4f})")
    print(f"End of exit parabola: ({Ex/ IN2M:.4f}, {Ey/ IN2M:.4f})")
    

    plt.plot(Nx/ IN2M, -Ny/ IN2M, "ro")
    plt.plot(Qx/ IN2M, -Qy/ IN2M, "go")
    plt.plot(Ex/ IN2M, -Ey/ IN2M, "bo")
    
    plt.plot(x_c, y_c, x_con_rad, y_con_rad, x_cone, y_cone, x_con, y_con, x_div, y_div, x_bell, y_bell)
    plt.title("Nozzle Contour")
    plt.xlabel("Axial Distance [in]")
    plt.ylabel("Radius [in]")
    plt.axis('equal')
    plt.show()



Dt = 2.3094013 #in
exp_ratio = 2.288
Lstar = 50
con_angle = 30
filename = 'chamber_contour'
theta_n = 20.88
theta_e = 14.6
Dc = 6
IN2M = 0.0254
Ac = ((((Dc * IN2M)/2)**2)*np.pi) #m^2
At= 0.027 #m^2
contract_ratio = Ac/At
print(f"Contraction Ratio: {contract_ratio:.4f}")


nozzle_contour(Dt, exp_ratio, Lstar, contract_ratio, con_angle, Dc, filename)
