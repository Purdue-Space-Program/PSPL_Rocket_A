import numpy as np
import matplotlib.pyplot as plt

def nozzle_contour(At, exp_ratio, Lstar, contract_ratio, con_angle, Dc, filename):
    # Unit conversions
    IN2M = 0.0254
    
    # Nozzle dimensions
    Lstar = Lstar * IN2M
    Dc = Dc * IN2M
    Rt = np.sqrt(At / np.pi)  # Throat radius
    Ae = At * exp_ratio # Exit area
    Re = np.sqrt(Ae / np.pi)  # Exit radius   
    percent_rao = 0.8
    Ln = percent_rao * ((np.sqrt(exp_ratio) - 1) * Rt) / np.tan(np.deg2rad(15)) # Nozzle length
    theta_n = np.deg2rad(-0.00001*exp_ratio**4 + 0.0017*exp_ratio**3 - 0.0775*exp_ratio**2 + 1.5769*exp_ratio + 16.805) # Initial bell angle
    theta_e = np.deg2rad(0.00001*exp_ratio**4 - 0.0015*exp_ratio**3 + 0.0665*exp_ratio**2 - 1.2849*exp_ratio + 17.895) # Final bell angle
    
    print(f"theta_n: {np.rad2deg(theta_n)}\n")
    if theta_e < 0:
        print("WHAT")
        theta_e = 0

    # Throat converging section arc
    con_angle = np.deg2rad(con_angle) # Convergence angle 
    theta_con = np.linspace(-np.pi/2 - con_angle, -np.pi/2)
    x_con = 1.5 * Rt * np.cos(theta_con)
    y_con = (1.5 * Rt * np.sin(theta_con)) + (1.5 * Rt) + Rt
    
    # Throat diverging section arc
    theta_div = np.linspace(-np.pi/2, theta_n - np.pi/2)
    x_div = 0.382 * Rt * np.cos(theta_div)
    y_div = (0.382 * Rt * np.sin(theta_div)) + (0.382 * Rt) + Rt
    x_div = np.delete(x_div, -1)
    y_div = np.delete(y_div, -1)
    
    # Bell section
    # pretty much all of this is on this: http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
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
    
    # this is a quadratic bezier curve where N, Q, and E are the 3 points being interpolated: https://en.wikipedia.org/wiki/B%C3%A9zier_curve#Quadratic_B%C3%A9zier_curves
    theta_bell = np.linspace(0, 1, num=20)
    x_bell = (((1 - theta_bell)**2) * Nx) + (2 * (1 - theta_bell) * theta_bell * Qx) + ((theta_bell**2) * Ex)
    y_bell = (((1 - theta_bell)**2) * Ny) + (2 * (1 - theta_bell) * theta_bell * Qy) + ((theta_bell**2) * Ey)

    # Chamber dimensions
    #Ac = contract_ratio * At
    #Rc = np.sqrt(Ac / np.pi)
    Ac = np.pi/4 * Dc**2
    Rc = Dc / 2
    contract_ratio = Ac / At

    # Converging cone
    Rcon = 1 * Rc
    x_length_cone = (Rc - Rcon * (np.sin(np.pi/2) - np.sin(np.pi/2 - con_angle)) - y_con[0]) / np.tan(-con_angle)
    x_cone = np.linspace(x_length_cone, 0) + x_con[0]
    y_cone = np.linspace(Rc - Rcon * (np.sin(np.pi/2) - np.sin(np.pi/2 - con_angle)) - y_con[0], 0) + y_con[0]
    x_cone = np.delete(x_cone, -1)
    y_cone = np.delete(y_cone, -1)

    # Converging radius
    theta_con_rad = np.linspace(np.pi/2, np.pi/2 - con_angle)
    x_con_rad = Rcon * np.cos(theta_con_rad) + x_cone[0] - Rcon * np.cos(np.pi/2 - con_angle)
    y_con_rad = Rcon * np.sin(theta_con_rad) - Rcon + Rc
    x_con_rad = np.delete(x_con_rad, -1)
    y_con_rad = np.delete(y_con_rad, -1)

    # Volume calculations
    V_con = np.pi * np.trapz(np.concatenate([y_con_rad, y_cone, y_con])**2, x=np.concatenate([x_con_rad, x_cone, x_con])) # Volume of converging section
    Vc = (Lstar * At) - V_con # Volume of cylindrical section 
    Lc = Vc / Ac # Length of chamber
    
    # Chamber
    x_c = np.linspace(-Lc, 0) + x_con_rad[0] 
    y_c = np.zeros(np.size(x_c)) + Rc
    x_c = np.delete(x_c, -1)
    y_c = np.delete(y_c, -1)

    x_arr = np.concatenate([x_c, x_con_rad, x_cone, x_con, x_div, x_bell])
    r_arr = np.concatenate([y_c, y_con_rad, y_cone, y_con, y_div, y_bell])

    chamber_length = x_arr[-1] - x_arr[0]
    chamber_diameter = 2 * np.sqrt(Ac / np.pi)

    print(f'Chamber diameter: {chamber_diameter / IN2M:.3f} in')
    print(f'Combustion chamber length: {Lc / IN2M:.3f} in')
    print(f'Total chamber length: {chamber_length / IN2M:.3f} in')
    print(f'Injector-to-throat length: {x_arr[0]*-1 / IN2M:.3f} in')
    print(f'Initial parabola angle: {np.rad2deg(theta_n):.2f} degrees')
    print(f'Exit parabola angle: {np.rad2deg(theta_e):.2f} degrees')
    
    z_arr = np.zeros(np.size(x_arr))
    nozzle = np.transpose(np.array([x_arr, r_arr, z_arr]))
    nozzle_in = nozzle / IN2M
    np.savetxt(filename, nozzle_in, delimiter=',')
    np.savetxt("nozzle_meters.csv", nozzle, delimiter=',')
    # np.savetxt("2d_meters", np.transpose(np.array([x_arr, r_arr])), delimiter=',')


    # Nozzle contour plot
    plt.plot(x_arr/IN2M, r_arr/IN2M, 'k', x_arr/IN2M, -r_arr/IN2M, 'k')
    
    print(f"Start: ({Nx:.4f}, {Ny:.4f})")
    print(f"Middle: ({Qx:.4f}, {Qy:.4f})")
    print(f"End: ({Ex:.4f}, {Ey:.4f})")
    
    # my attempt to figure out whats wrong with this: https://www.desmos.com/calculator/fdnhc6le8s
    # When you plot this and zoom in you'll see that these are not in the order
    # red, green, blue from left to right even through they should be.
    # My guess is the middle point is wrong.
    plt.plot(Nx, Ny, "ro")
    plt.plot(Qx, Qy, "go")
    plt.plot(Ex, Ey, "bo")
    
    plt.plot(x_c, y_c, x_con_rad, y_con_rad, x_cone, y_cone, x_con, y_con, x_div, y_div, x_bell, y_bell)
    plt.title("Nozzle Contour")
    plt.xlabel("Axial Distance [in]")
    plt.ylabel("Radius [in]")
    plt.axis('equal')
    plt.show()
    


    return chamber_diameter, Lc, chamber_length, x_arr[0]*-1

At = 0.00270322
exp_ratio = 3
Lstar = 15
contract_ratio = 2.288
con_angle = 30
filename = 'Pathfinder Rocket Nozzle Contour3'
theta_n = 20.88
theta_e = 14.6
Dc = 3.494



nozzle_contour(At, exp_ratio, Lstar, contract_ratio, con_angle, Dc, filename)
