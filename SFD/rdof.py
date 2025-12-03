import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator
from scipy.integrate import solve_ivp
from SFD import getRocketSections, returnMassModel, returnFinCP, returnFinSD
from ambiance import Atmosphere

plt.rcParams['lines.linewidth'] = 0.5
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Libre Franklin', 'Arial']
plt.rcParams['axes3d.mouserotationstyle'] = 'azel'
plt.rcParams["axes.grid"] = True
plt.rcParams['legend.fontsize'] = 'small'

PRIMARY_LINE_COLOR = "#DDB945"
SECONDARY_LINE_COLOR = "#8E6F3E"

# Conversion
LB_TO_KG = 0.453592
KG_TO_LB = 1 / LB_TO_KG

FT_TO_M = 0.3048
M_TO_FT = 1 / FT_TO_M

IN_TO_M = 0.0254
M_TO_IN = 1 / IN_TO_M

N_TO_LBF = 0.224809
NEWTON_METER_TO_FOOT_POUND = 0.737562 # 1 N*m = 0.737562 ft*lbf


# Timestep parameters
dt = 0.01
SIM_TIME = 600 # sec
t_vec = np.arange(0, SIM_TIME, dt)

# Rocket parameters
M = 180.98 * LB_TO_KG # kg, subtracting nosecone mass

OD = 8.625 * IN_TO_M  # m
radius = OD / 2  # m

phi_init = 0 # degrees
L = 245* IN_TO_M - 3.59 * FT_TO_M # rocket length (m) subtract nosecone length
CGpos = 2.45  # CG position from aft (m), adjusted for no nosecone
CGtoTip = L - CGpos  # CG position from tip (m)
CGtoCenter = L/2 - CGpos  # CG position from center (m)
I = 1/4 * M * pow(radius, 2) + 1/12 * M * pow(L, 2) + M*pow(CGtoCenter,2) # kg*m^2

rocketLongArea = L * OD # m^2, area of rectangle
rocketShortArea = np.pi * pow(radius, 2)  # m^2

CD_longSide = 2  # drag coefficient of the cylinder perpendicular to flow
CD_shortSide = 0.82  # drag coefficient of the cylinder parallel to flow
CD_cylinder = 1.1
C_rotational = 20

CD_drogue = .97 # drag coefficient of drogue parachute
CD_main = 2.2 # drag coefficient of main parachute

# Initial conditions
APOGEE = 40_000 # ft
WORST_CASE_ALT = 10_000 # ft
worst_case_alt = APOGEE * FT_TO_M
vx_init = 131.234 * FT_TO_M
t0 = [0, 0, vx_init, 0, worst_case_alt]  # [phi, omega, vx, vy, altitude]

# Create chute area piecewise function
drogue_radius = 4 * FT_TO_M # m
drogue_area = np.pi * pow(drogue_radius, 2)  # m^2
inflation_time_n = .2 # sec
drogue_opening_time = inflation_time_n * (drogue_radius*2) / vx_init # seconds

main_radius = 9 * FT_TO_M # m
main_area = np.pi * pow(main_radius, 2)  # m^2
main_opening_time = 2 # seconds

main_altitude_time = 0 # set this variable to the time at which the main chute opens
MAIN_DEPLOYMENT_ALTITUDE = 1000 # ft

# Fin stability derivative and center of pressure from SFD
FIN_SD = returnFinSD(OD)
FIN_CP = returnFinCP(L-0.5, L)

def chute_drag_area(t, alt):
        if alt < MAIN_DEPLOYMENT_ALTITUDE * FT_TO_M: 
            return main_area * CD_main
        return drogue_area * CD_drogue

# Helper functions

def magnitude(vec):
    return np.sqrt(vec[0]**2 + vec[1]**2)

def sin(theta): return np.sin(theta)
def cos(theta): return np.cos(theta)
def atan(y, x): return np.arctan2(y, x)

def normalize(vec):
    mag = magnitude(vec)
    if mag == 0: return vec
    return vec / mag

def rho(alt): return np.exp(
        4.88158e-18 * alt**4
        - 1.808e-13 * alt**3
        + 2.432e-11 * alt**2
        - 9.693e-5 * alt
        + 0.1922
    ) # air density at a certain altitude 
# def rho(alt): return Atmosphere(alt).density[0]  # kg/m^3, air density at a certain altitude

def freeStreamAngle(vel):
    v_g1 = vel[0]
    v_g2 = vel[1]
    angle = atan(v_g2, v_g1)
    return angle

def transform_gTor(g_vec, phi):
    g1, g2 = g_vec
    r1 = g1 * cos(phi) + g2 * sin(phi)
    r2 = -g1 * sin(phi) + g2 * cos(phi)
    return np.array([r1, r2])

def dragArea(vx, vy, phi): # Rocket body drag area, dot product of velocity & normal dir
    normal_dir = [-sin(phi), cos(phi)]  # normal direction to rocket axis
    normalized_v = np.array([vx, vy]) / magnitude([vx, vy])  # normalized velocity vector
    cylinder_side_area = L * OD
    return np.dot(normalized_v, normal_dir) * cylinder_side_area * CD_cylinder

# Force calculations
def F_g():
    return np.array([0, -M * 9.81])  # N, gravitational force

def F_s(t, alt, phi, omega, vx, vy):
    v = [vx, vy]
    opening_factor = (1 - t) if t <= 1 else 0
    windSpeedMax = 80 # m/s
    v[0] -= windSpeedMax * opening_factor
    rocketVelMag = CGtoTip * omega # m/s, velocity of the rocket tip
    chuteVel = np.array(v + rocketVelMag * np.array([-sin(phi), cos(phi)]))  # m/s, velocity of the rocket tip
    magVel = magnitude(chuteVel)  # m/s, magnitude of the velocity vector
    V_norm = normalize(chuteVel)  # normalized velocity vector

    chuteDragArea = chute_drag_area(t, alt)  # m^2
    Fs_mag = 0.5 * rho(alt) * (magVel)**2 * chuteDragArea
    Fs_mag *= 1 + opening_factor * 0.8
    return Fs_mag * -V_norm # N, parachute drag force

def F_bodydrag(alt, vx, vy, omega, phi):

    v_norm = normalize([vx, vy])  # normalized velocity vector
    drag_area = dragArea(vx, vy, phi)  # m^2
    magVel = magnitude([vx, vy])  # m/s, magnitude of the velocity vector

    Fbd_mag = 0.5 * rho(alt) * magVel**2 * drag_area  # N, drag forceorce
    return Fbd_mag * -v_norm

def M_bodydrag(alt, vx, vy, phi, dist):
    velAngle = freeStreamAngle([vx, vy])  # degrees

    longSide_area = CGtoTip * OD
    negativeSide_area = CGtoCenter * OD
    drag_area = longSide_area - negativeSide_area  # m^2

    magVel = magnitude([vx, vy])
    moment_angle = phi - velAngle
    moment = 0.5 * rho(alt) * magVel**2 * drag_area * sin(moment_angle) * dist  # N*m, moment due to drag force
    return moment  # N*m, moment due to body drag

def M_chute(t, alt, omega, vx, vy, phi, dist):
    # return moment  # N*m, moment due to parachute drag force
    F_s_shear = transform_gTor(F_s(t, alt, phi, omega, vx, vy), phi)[1]
    return F_s_shear * dist

def rotationalDrag(alt, omega, tipDistance, centerDistance):
    coefficient = (CD_cylinder * rho(alt) * OD * omega**2)/8
    distances = (tipDistance**4 + centerDistance**4)
    moment = -np.sign(omega) * coefficient * distances  # N*m, moment due to rotational drag
    return moment

def calcMoments(t, alt, omega, vx, vy, phi, pos):
    # Pos is position from aft in m
    posToCenter = L/2 - pos
    posToTip = L - pos

    Mbd = M_bodydrag(alt, vx, vy, phi, posToCenter)  # N*m, moment due to body drag
    Mrot = rotationalDrag(alt, omega, posToTip, posToCenter)  # N*m, moment due to rotational drag
    Mchute = M_chute(t, alt, omega, vx, vy, phi, posToTip)  # N*m, moment due to parachute drag

    return Mbd + Mrot + Mchute

def calc_alpha(t, alt, vx, vy, phi, omega):
    M_sum = calcMoments(t, alt, omega, vx, vy, phi, CGpos)

    alpha = M_sum / I  # rad/s^2, angular acceleration
    return alpha 

def force_system(t, state):
    # unpack state vector
    phi, omega, vx, vy, alt = state

    # Calculate forces
    Fg = F_g()
    Fs = F_s(t, alt, phi, omega, vx, vy)  # N, force due to parachute
    Fbd = F_bodydrag(alt, vx, vy, omega, phi)

    F_total = Fg + Fs + Fbd

    acceleration = F_total / M  # m/s^2
    ax, ay = acceleration

    # Calculate torque
    alpha = calc_alpha(t, alt, vx, vy, phi, omega)  # rad/s^2

    return [omega, alpha, ax, ay, vy]

def get_axial_compression(omega, mass_model,  F_total):

    axial_total_forces = -F_total[0]
    axial_accceleration = axial_total_forces/M

    sec_lengths = np.linspace(0, L, len(mass_model))
    CGtoSections = sec_lengths - CGpos
    centripetal_force = np.abs(mass_model * pow(omega, 2) * CGtoSections)

    # axial_forces = np.abs(axial_total_forces) - np.abs(mass_model * axial_accceleration)
    axial_forces = axial_total_forces + centripetal_force - mass_model * axial_accceleration  # N, axial forces across rocket body
    return axial_forces

def get_shear_forces(phi, alpha, mass_model, F_chute, F_bodydrag, F_total):

    sec_lengths = np.linspace(0, L, len(mass_model))
    angular_acceleration_linearized = alpha * (sec_lengths - CGpos)
    accel_shear = F_total[1]/M
    mass_model_slices = np.diff(mass_model, prepend=0)  # Mass of each element, not cumulative
    sec_Fbds = (F_bodydrag[1]/L) * sec_lengths # Distributed body drag over length of rocket

    inertial_forces = np.cumsum(mass_model_slices * (accel_shear + angular_acceleration_linearized))
    shear_forces = -inertial_forces
    F_gs = mass_model * -9.81 # N, gravitational force on each segment
    F_gs_r = [transform_gTor([0, Fg], phi) for Fg in F_gs]  # N, gravitational force in rocket frame
    F_g_shears = [g[1] for g in F_gs_r]  # N, shear force due to gravity in rocket frame
    shear_forces += sec_Fbds + F_g_shears
    shear_forces[-12:] += F_chute[1]

    return shear_forces

def get_bending_moments(dy, shear_forces):

    return np.cumsum(shear_forces) * dy # N*m, bending moments across rocket body

def plot(subplot, x_vec, y_vecs, title, x_label, y_label, legends=['']):
    BACKGROUND_COLOR = "#e8e8e8"
    LINE_COLORS = ["#daaa00", "#8E6F3E", "#DDB945", "#CFB991"]

    subplot.set_facecolor(BACKGROUND_COLOR)
    if len(y_vecs) == len(x_vec): # Single line
        y_vecs = [y_vecs]
    for i, y_vec in enumerate(y_vecs):
        subplot.plot(x_vec, y_vec, color=LINE_COLORS[i], label=legends[i], linewidth=1)
    subplot.set_title(title)
    subplot.set_xlabel(x_label)
    subplot.set_ylabel(y_label)
    # subplot.grid()

    if legends != ['']:
        subplot.legend()

def default():
    global t_vec, dy, mass_model

    # Use solve_ivp to solve ODE system
    t_span = (0, SIM_TIME)
    solution = solve_ivp(force_system, t_span, t0, t_eval=t_vec, method='RK45')
    
    # Unpack solution
    phi = solution.y[0]
    omega = solution.y[1]
    vx = solution.y[2]
    vy = solution.y[3]
    altitude = solution.y[4]

    # If t_vec longer than solution length, match lengths
    if len(t_vec) > len(solution.t):
        t_vec = t_vec[:len(solution.t)]

    # Calculate forces across t_vec
    Fg = np.array([F_g() for _ in t_vec]) * N_TO_LBF
    Fbd = np.array([F_bodydrag(altitude[i], vx[i], vy[i], omega[i], phi[i]) for i in range(len(vx))]) * N_TO_LBF  # lbf, body drag force
    Fs = np.array([F_s(t_vec[i], altitude[i], phi[i], omega[i], vx[i], vy[i]) for i in range(len(vx))]) * N_TO_LBF  # lbf, parachute drag force
    Fg_magnitude = np.array([magnitude(f) for f in Fg])
    Fbd_magnitude = np.array([magnitude(f) for f in Fbd])
    Fs_magnitude = np.array([magnitude(f) for f in Fs])

    # Sum forces and convert to Rocket frame
    F_total = Fg + Fbd + Fs # N, total force
    F_total_r = np.array([transform_gTor(F_total[i], phi[i]) for i in range(len(F_total))])  # lbf, total force in rocket frame

    # Calculate moments
    alphas = np.array([calc_alpha(t_vec[i], altitude[i], vx[i], vy[i], phi[i], omega[i]) for i in range(len(vx))])  # rad/s^2, angular acceleration
    M_total = alphas * I


    dy = 5/1000 # m
    mass_model = returnMassModel(dy)
    mass_model = np.cumsum(mass_model)  # kg*m^2, mass model
    length_vec = np.linspace(0, L, len(mass_model))  # m, position vector from aft to tip

    F_chute_r = np.array([transform_gTor(Fs[i], phi[i]) for i in range(len(Fs))]) / N_TO_LBF  # N, parachute drag force in rocket frame
    F_bodydrag_r = np.array([transform_gTor(Fbd[i], phi[i]) for i in range(len(Fbd))]) / N_TO_LBF 

    import time
    start_time = time.time()

    # make two subplots
    shear_bending_fig, (axial3d, shear3d, bending3d) = plt.subplots(1, 3, subplot_kw={'projection': '3d'})
    shear_bending_fig.suptitle("Internal Loads over Rocket Length over Time", fontsize=16)

    #3D arrays for shear and bending
    shear_forces_over_time = []
    bending_moments_over_time = []
    axial_forces_over_time = []     

    shear_bending_dt = .1
    shear_bending_start = 0
    shear_bending_end = 5
    downsample_factor = 5
    t_vec_shear_bending = np.arange(shear_bending_start, shear_bending_end, shear_bending_dt)  # seconds

    for t in t_vec_shear_bending: # iterate over first __ seconds of flight
        mi = int(t/dt) #
        
        shear = get_shear_forces(phi[mi], alphas[mi], 
                                mass_model, F_chute_r[mi], F_bodydrag_r[mi], F_total_r[mi]/N_TO_LBF)
        shear_forces_over_time.append(shear)  # N, shear forces across rocket body
    
        # M0 = M_total[mi] # N*m, moment at CG position
        bending = get_bending_moments(dy, shear)
        # if AoA[mi] > np.pi: bending *= -1
        bending_moments_over_time.append(bending) 

        axial_forces = get_axial_compression(omega[mi], mass_model, F_total_r[mi]/N_TO_LBF)
        axial_forces_over_time.append(axial_forces)  # N, axial forces across rocket body
   
    # plot surface on 3d plot
    length_3d, time_3d = np.meshgrid(np.array(length_vec) * M_TO_FT, np.array(t_vec_shear_bending), indexing='ij')

    bending3d.plot_surface(time_3d, length_3d, np.array(bending_moments_over_time).T*NEWTON_METER_TO_FOOT_POUND, cmap='viridis', rstride=downsample_factor, cstride=downsample_factor, alpha=0.9, edgecolor='none')
    bending3d.set_xlabel("Time (s)")
    bending3d.set_ylabel("Length from Aft (ft)")
    bending3d.set_zlabel("Bending Moment (ft-lbs)")
    bending3d.yaxis.set_major_locator(MultipleLocator(3))
    bending3d.set_title("Bending Moment")
    bending3d.invert_xaxis() 


    shear3d.plot_surface(time_3d, length_3d, np.array(shear_forces_over_time).T*N_TO_LBF, cmap='viridis', rstride=downsample_factor, cstride=downsample_factor, alpha=0.9, edgecolor='none')
    shear3d.set_xlabel("Time (s)")
    shear3d.set_ylabel("Length from Aft (ft)")
    shear3d.set_zlabel("Shear Force (lbf)")
    shear3d.yaxis.set_major_locator(MultipleLocator(3))
    shear3d.set_title("Shear Force")
    shear3d.invert_xaxis() 

    axial3d.plot_surface(time_3d, length_3d, np.array(axial_forces_over_time).T * N_TO_LBF, cmap='plasma', rstride=downsample_factor, cstride=downsample_factor, alpha=0.9, edgecolor='none')
    axial3d.set_xlabel("Time (s)")
    axial3d.set_ylabel("Length from Aft (ft)")
    axial3d.set_zlabel("Axial Compression (lbf)")
    axial3d.set_title("Axial Compression")
    axial3d.yaxis.set_major_locator(MultipleLocator(3))
    axial3d.invert_xaxis() 

    # Find max. axial compression, shear force & bending moment at each rocket point
    worst_axial = np.max(np.abs(np.array(axial_forces_over_time)), axis=0)
    worst_shear = np.max(np.abs(np.array(shear_forces_over_time)), axis=0)
    worst_bending = np.max(np.abs(np.array(bending_moments_over_time)), axis=0)
    # Find worst shock load
    max_chute_load_time = np.argmax(Fs_magnitude)
    max_chute_axial_load = transform_gTor(Fs[max_chute_load_time], phi[max_chute_load_time])

    end_time = time.time()
    print(f"3D Graph Computation time: {end_time - start_time:.2f} seconds")
    # Save internal loads to excel
    import output_formatter
    import pandas as pd
    
    input_xlsx = pd.ExcelFile('inputs.xlsx', engine='openpyxl')
    rocket_dict = getRocketSections(input_xlsx)
    output_formatter.outputFormatter(shear=worst_shear, bending=worst_bending, axial=worst_axial, shock_load=max_chute_axial_load, mass=M, length=L, rocket_dict=rocket_dict)

    plt.show() # Show 3D plots

    # Create figure with 3 subplots: Axial, shear & bending
    internal_fig, (axial_plot, shear_plot, bending_plot) = plt.subplots(1, 3, figsize=(12, 4))
    internal_fig.suptitle("Worst-Case Internal Rocket Loads during Recovery", fontsize=16)
    WINDOW_COLOR = "white"
    internal_fig.patch.set_facecolor(WINDOW_COLOR)

    plot(axial_plot, length_vec * M_TO_FT, worst_axial * N_TO_LBF, 'Axial Compression', 'Length from Aft (ft)', 'Axial Compression (lbf)')
    plot(shear_plot, length_vec * M_TO_FT, worst_shear * N_TO_LBF, 'Shear Force', 'Length from Aft (ft)', 'Shear Force (lbf)')
    plot(bending_plot, length_vec * M_TO_FT, worst_bending * NEWTON_METER_TO_FOOT_POUND, 'Bending Moment', 'Length from Aft (ft)', 'Bending Moment (ft-lbs)')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

    altitude *= M_TO_FT
    vx *= M_TO_FT
    vy *= M_TO_FT

    # Plot a 2x4 grid of subplots:
    fig, axs = plt.subplots(2, 4, figsize=(12, 8))
    fig.suptitle(f"RDOF - {drogue_radius*M_TO_FT:.1f} ft. drogue, {main_radius*M_TO_FT:.1f} ft. main", fontsize=16)
    fig.patch.set_facecolor(WINDOW_COLOR)

    # Subplot 1: Velocity x and y vs time
    plot(axs[0,0], t_vec, [vx, vy], 'Velocity vs Time', 'Time (s)', 'Veocity (ft/s)', legends=['Vx', 'Vy'])
   
    # Subplot 2: Individual forces vs time
    plot(axs[0,1], t_vec, [Fg_magnitude, Fbd_magnitude, Fs_magnitude], 'Global Forces vs Time', 'Time (s)', 'Force (lbf)', legends=['Weight', 'Body Drag', 'Chute Drag'])

    # Subplot 3: Phi and velocity angle vs time
    phi_vertical = np.pi/2 - phi
    plot(axs[1,0], t_vec, [phi_vertical * 180/
                           np.pi], 'Rocket Angle vs Time', 'Time (s)', 'Angle from Vertical (degrees)', legends=['Rocket Angle'])

    # Subplot 4: Axial Compression across rocket body
    plot(axs[1,1], length_vec * M_TO_FT, worst_axial * N_TO_LBF, 'Axial Compression', 'Length from Aft (ft)', 'Axial Compression (lbf)')
    # Subplot 5: Rocket frame forces vs time
    plot(axs[0,2], t_vec, [F_total_r[:, 0], F_total_r[:, 1]], 'Rocket Frame Forces vs Time', 'Time (s)', 'Force (lbf)', legends=['Axial', 'Shear'])
    # Subplot 6: Internal bending moments
    plot(axs[1,2], length_vec * M_TO_FT, worst_bending * NEWTON_METER_TO_FOOT_POUND, 'Bending Moment', 'Length from Aft (ft)', 'Bending Moment (ft-lbs)')
    # Subplot 7: Altitude vs time
    plot(axs[0,3], t_vec, altitude, 'Altitude vs Time', 'Time (s)', 'Altitude (ft)')

    # Animation
    from matplotlib.animation import FuncAnimation

    # Subplot 8: Rocket orientation animation
    animSubplot = plt.subplot(2, 4, 8, polar=True)

    title = animSubplot.text(0.5,0.1, "Rocket Orientation at t=0.0s",
                    transform=animSubplot.transAxes, ha="center", bbox={'facecolor': 'w', 'alpha': 0.75, 'pad': 2})

    polarAxis = animSubplot
    # remove numbers from the polar axis
    polarAxis.set_xticklabels([])
    polarAxis.set_yticklabels([])
    polarAxis.set_xlabel(None)
    polarAxis.set_ylabel(None)

    # Make 0 degrees at the top
    polarAxis.set_theta_zero_location('N')
    # Make np.pi/2 degrees at the right
    polarAxis.set_theta_direction(-1)

    # Draw line representing rocket orientation
    line, = polarAxis.plot([np.pi/2, np.pi/2], [0, CGtoTip], lw=2, color='black')
    bottomLine, = polarAxis.plot([1.5*np.pi, 1.5 * np.pi], [0, CGpos], lw=2, color='black')

    anim_speed_factor = 3

    def update(i):

        i = int(i * anim_speed_factor)

        string = f"t = {t_vec[i]:.1f}s\nAltitude = {altitude[i]:.1f}ft"
        title.set_text(string)

        rocket_angle = phi_vertical[i] # degrees
        rocket_angle_rad = rocket_angle

        # update line to point in rocket_angle_rad direction, and length of line is CGtoTip

        line.set_data([rocket_angle_rad, rocket_angle_rad], [0, CGtoTip])
        bottomLine.set_data([rocket_angle_rad+np.pi, rocket_angle_rad+np.pi], [0, CGpos])

        return line, bottomLine, title

    anim = FuncAnimation(fig, update, frames=len(t_vec)//anim_speed_factor, blit=True, interval=1, repeat=True)

    # Save the figure as a PNG file
    # fig.savefig('RDOF.png', dpi=300, bbox_inches='tight')

    # Find ground impact speed
    groundIndex = np.where(altitude <= 0)
    if groundIndex[0].size > 0:
        ind = groundIndex[0][0]
        v = magnitude([vx[ind], vy[ind]])
        print(f"Rocket velocity at ground: {v:.2f} ft/s")
    else:
        print("Rocket doesn't reach ground given timeframe/conditions.")

    # Show animation
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

def getDataForChuteSize(drogue, main): # drogue and main chute radius (m)
    global drogue_area, main_area, main_altitude_time

    main_altitude_time = 0

    drogue_area = np.pi * pow(drogue, 2)
    main_area = np.pi * pow(main, 2)

    # Use solve_ivp to solve ODE system
    t_span = (0, SIM_TIME)
    solution = solve_ivp(force_system, t_span, t0, t_eval=t_vec, method='RK45')
    
    # Unpack solution
    phi = solution.y[0]
    omega = solution.y[1]
    vx = solution.y[2]
    vy = solution.y[3]
    altitude = solution.y[4]

    Fg = np.array([F_g() for _ in t_vec])
    Fbd = np.array([F_bodydrag(altitude[i], vx[i], vy[i], omega[i], phi[i]) for i in range(len(vx))])  # lbf, body drag force
    Fs = np.array([F_s(t_vec[i], altitude[i], phi[i], omega[i], vx[i], vy[i]) for i in range(len(vx))]) # lbf, parachute drag force
 
    # Sum forces and convert to Rocket frame
    F_total = Fg + Fbd + Fs  # N, total force
    F_total_r = np.array([transform_gTor(F_total[i], phi[i]) for i in range(len(F_total))])  # lbf, total force in rocket frame

    groundIndex = np.where(altitude <= 0)[0]
    mainIndex = np.where(altitude <= MAIN_DEPLOYMENT_ALTITUDE * FT_TO_M)[0]
    vImpact = tImpact = tChute = -1
    if groundIndex.size > 0:
        ind = groundIndex[0]
        vImpact = magnitude([vx[ind], vy[ind]])
        tImpact = t_vec[ind]

    if mainIndex.size > 0:
        tChute = t_vec[mainIndex[0]]

    return F_total_r * N_TO_LBF, vImpact * M_TO_FT, tImpact, tChute

def forcesVsChuteSize(drogues, mains): # chute sizes IN FEET!!
    fig, axs = plt.subplots(nrows=len(drogues), ncols=len(mains), sharey='col', squeeze=False)
    fig.suptitle('Chute Diameters vs. Rocket Forces & Impact Velocity', fontsize=16)

    for i, drogue in enumerate(drogues): # rows
        for j, main in enumerate(mains): # columns
            F_r, v_impact, t_impact, t_chute = getDataForChuteSize(drogue * FT_TO_M, main * FT_TO_M)
            plot = axs[i, j]
            plot.plot(t_vec, F_r[:, 0], label='Axial', color=PRIMARY_LINE_COLOR)
            plot.plot(t_vec, F_r[:, 1], label='Shear', color=SECONDARY_LINE_COLOR)
            
            # Annotate chute deployment & impact timestamps
            if t_chute != -1: 
                plot.annotate("Main", xy=(t_chute, 0), horizontalalignment='right', xytext=(t_chute, 5))
                plot.plot([t_chute], [0], marker='.', color='red', linestyle='None')
            if t_impact != -1: 
                plot.plot([t_impact], [0], marker='.', color='red', linestyle='None')
                plot.annotate("Impact", xy=(t_impact, 0), xytext=(t_impact, 5))
            
            if v_impact >= 0 :plot.set_title(f"{drogue*2:.0f} ft. drogue, {main*2:.0f} ft. main - Vimp: {v_impact:.1f} ft/s")
            else: plot.set_title(f"{drogue*2:.0f} ft. drogue, {main*2:.0f} ft. main - Vimp: N/A")
            if i == len(drogues) - 1:  plot.set_xlabel('Time (s)')
            else: plot.set_xticklabels([])
            if j == 0: plot.set_ylabel('Force (lbf)')
            plot.legend(fontsize="small")
            plot.grid()
    
    # plt.show()


default()

drogues = np.array([3, 4,5])
mains = np.array([6, 8, 10])
# forcesVsChuteSize(drogues, mains)