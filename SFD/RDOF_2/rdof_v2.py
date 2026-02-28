import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.ticker import MaxNLocator, MultipleLocator
from scipy.integrate import solve_ivp
from scipy.io import savemat

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


import vehicle_parameters
import vehicle_parameters_functions
from vehicle_parameters import parameters, rocket_dict_dry, rocket_dict_wet
import vehicle_main
import constants as c

try:
    # from SFD import parseWind
    # from SFD import sfd
    from SFD import loads
except ModuleNotFoundError:
    # import parseWind
    # import sfd
    import loads




# Mass Model
def get_mass_model(rocket_dict, parachute_mass):
    '''
    rocket_dict: Dictionary of rocket components
    parachute_mass: Mass of parachute [kg]
    linear_density_array: Array of linear density along the rocket [kg / m]
    length_along_rocket_linspace: Array of length along rocket [m]
    '''
    num_points = 500
    length_along_rocket_linspace = np.linspace(rocket_dict["engine"]["bottom_distance_from_aft"], rocket_dict["nosecone"]["bottom_distance_from_aft"] + rocket_dict["nosecone"]["length"], num_points)
    dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]  # [m]

    linear_density_array = np.zeros(num_points)

    for component in rocket_dict:
        mass = rocket_dict[component]["mass"]
        if component == "recovery_bay":
            mass -= parachute_mass
        bottom_distance_from_aft = rocket_dict[component]["bottom_distance_from_aft"]
        length = rocket_dict[component]["length"]
        linear_density = mass / length
        for index, length_along_rocket in enumerate(length_along_rocket_linspace):
            above_component_bottom = length_along_rocket >= bottom_distance_from_aft
            below_component_top = length_along_rocket <= (bottom_distance_from_aft + length)
            
            if (above_component_bottom and below_component_top):
                linear_density_array[index] += linear_density
    recovery_bay_start = rocket_dict["recovery_bay"]["bottom_distance_from_aft"]
    index = int(recovery_bay_start / dx)
    linear_density_array = linear_density_array[:index]
    length_along_rocket_linspace = length_along_rocket_linspace[:index]

    return linear_density_array, length_along_rocket_linspace







############ FIX main_altitude_time #################








def Calculate_Chute_Drag_Area(time, altitude, MAIN_DEPLOYMENT_ALTITUDE, drogue_drag_area, drogue_opening_time, main_drag_area, main_opening_time, main_altitude_time):
    if altitude < MAIN_DEPLOYMENT_ALTITUDE * c.FT2M:
        if main_altitude_time is None:
            main_altitude_time = time
            chute_drag_area = drogue_drag_area
        elif time - main_altitude_time < main_opening_time:
            opening_factor = (1 - (time-main_altitude_time)/main_opening_time)
            chute_drag_area = max(drogue_drag_area, main_drag_area * ((time - main_altitude_time) / main_opening_time)) * (1+opening_factor)
        else:
            chute_drag_area = main_drag_area
    else:
        if time < drogue_opening_time:
            opening_factor = (1-time/drogue_opening_time)
            chute_drag_area = drogue_drag_area * (time / drogue_opening_time) * (1+opening_factor)
        else:
            chute_drag_area = drogue_drag_area

    return(chute_drag_area, main_altitude_time)

# Helper functions

def magnitude(vector):
    return np.sqrt(vector[0]**2 + vector[1]**2)

def sin(theta): return np.sin(theta)
def cos(theta): return np.cos(theta)
def atan(y, x): return np.arctan2(y, x)

def normalize(vec):
    mag = magnitude(vec)
    if mag == 0: return vec
    return vec / mag

# Air density at a certain altitude
def Get_Air_Density(altitude): return np.exp(
        4.88158e-18 * altitude**4
        - 1.808e-13 * altitude**3
        + 2.432e-11 * altitude**2
        - 9.693e-5 * altitude
        + 0.1922
    ) # air density at a certain altitude

# Transform from global reference frame to rocket reference frame
def transform_gTor(g_vec, phi):
    g1, g2 = g_vec
    r1 = g1 * cos(phi) + g2 * sin(phi)
    r2 = -g1 * sin(phi) + g2 * cos(phi)
    return np.array([r1, r2])

# Rocket body drag area, dot product of velocity & normal direction
def dragArea(vx, vy, phi, total_length, total_mass, rocket_outer_diameter, CD_cylinder):
    normal_dir = [-sin(phi), cos(phi)]  # normal direction to rocket axis
    normalized_v = np.array([vx, vy]) / magnitude([vx, vy])  # normalized velocity vector
    cylinder_side_area = total_length * rocket_outer_diameter
    return np.abs(np.dot(normalized_v, normal_dir) * cylinder_side_area * CD_cylinder)

# Force calculations
def Calculate_Force_Gravity(total_mass):
    return np.array([0, -total_mass * 9.81])  # N, gravitational force

# Parachute drag force
def Calculate_Force_Parachute_Drag(time, altitude, phi, omega, velocity_x, velocity_y, CGtoTip, MAIN_DEPLOYMENT_ALTITUDE, drogue_drag_area, drogue_opening_time, main_drag_area, main_opening_time, main_altitude_time):
    velocity = [velocity_x, velocity_y] # Velocity vector
    gust_factor = (1 - time) if time <= 1 else 0
    windSpeedMax = loads.parseWind.percentile_90_wind_gust_speed # 12 # m/s, disabled this for now
    velocity[0] -= windSpeedMax * gust_factor # Adding wind gusts in x direction
    rocketVelMag = CGtoTip * omega # m/s, velocity of the rocket tip, arc length equation / time (radius * angle / time)
    chuteVel = np.array(velocity + rocketVelMag * np.array([-sin(phi), cos(phi)]))  # m/s, velocity of the rocket tip
    velocity_magnitude = magnitude(chuteVel)  # m/s, magnitude of the parachute velocity vector
    V_norm = normalize(chuteVel)  # normalized parachute velocity vector

    chute_drag_area, main_altitude_time = Calculate_Chute_Drag_Area(time, altitude, MAIN_DEPLOYMENT_ALTITUDE, drogue_drag_area, drogue_opening_time, main_drag_area, main_opening_time, main_altitude_time)  # m^2
    Fs_mag = 0.5 * Get_Air_Density(altitude) * (velocity_magnitude)**2 * chute_drag_area
    # Fs_mag *= 1 + opening_factor * 0.8
    force_parachute_drag = Fs_mag * -V_norm
    return (force_parachute_drag, main_altitude_time) # N, parachute drag force

# Body drag force
def Calculate_Force_Body_Drag(alt, velocity_x, velocity_y, phi, total_length, total_mass, rocket_outer_diameter, CD_cylinder):
    v_normalized = normalize([velocity_x, velocity_y])  # normalized velocity vector
    drag_area = dragArea(velocity_x, velocity_y, phi, total_length, total_mass, rocket_outer_diameter, CD_cylinder)  # m^2
    magVel = magnitude([velocity_x, velocity_y])  # m/s, magnitude of the velocity vector

    Fbd_mag = 0.5 * Get_Air_Density(alt) * magVel**2 * drag_area  # N, drag force magnitude
    return Fbd_mag * -v_normalized

# Calculate moment due to body drag of rocket
def Calculate_M_body_drag(alt, velocity_x, velocity_y, phi, distance, total_length, total_mass, rocket_outer_diameter, CD_cylinder):
    F_body_drag = Calculate_Force_Body_Drag(alt, velocity_x, velocity_y, phi, total_length, total_mass, rocket_outer_diameter, CD_cylinder) # Body drag force
    F_body_drag_shear = transform_gTor(F_body_drag, phi)[1] # Shear component of body drag force in rocket reference frame
    moment_due_to_body_drag = F_body_drag_shear * distance 
    return moment_due_to_body_drag

# Calculate moment due to parachute drag of rocket
def Calculate_M_chute(time, altitude, omega, vx, vy, phi, dist, CGtoTip, MAIN_DEPLOYMENT_ALTITUDE, drogue_drag_area, drogue_opening_time, main_drag_area, main_opening_time, main_altitude_time):
    F_s_shear = transform_gTor(Calculate_Force_Parachute_Drag(time, altitude, phi, omega, vx, vy, CGtoTip, MAIN_DEPLOYMENT_ALTITUDE, drogue_drag_area, drogue_opening_time, main_drag_area, main_opening_time, main_altitude_time)[0], phi)[1] # Shear component of parachute
    return F_s_shear * dist

# Calculate rotational drag, it is a moment
def calculate_rotational_drag(alt, omega, total_length, CD_cylinder, rocket_radius):
    moment = -(Get_Air_Density(alt)* CD_cylinder * rocket_radius * pow(total_length,4)*omega*np.abs(omega))/12
    return moment

# Calculate total moment on rocket at a given time, altitude, velocity, orientation, and position along the rocket
def calculate_moments(time, altitude, omega, vx, vy, phi, pos, CD_cylinder, rocket_radius, rocket_outer_diameter, total_length, total_mass, CGtoTip, MAIN_DEPLOYMENT_ALTITUDE, drogue_drag_area, drogue_opening_time, main_drag_area, main_opening_time, main_altitude_time):
    # Pos is position from aft in m
    position_to_center = total_length/2 - pos
    posToTip = total_length - pos

    M_body = Calculate_M_body_drag(altitude, vx, vy, phi, position_to_center, total_length, total_mass, rocket_outer_diameter, CD_cylinder)  # N*m, moment due to body drag
    M_rot = calculate_rotational_drag(altitude, omega, total_length, CD_cylinder, rocket_radius)  # N*m, moment due to rotational drag
    M_chute = Calculate_M_chute(time, altitude, omega, vx, vy, phi, posToTip, CGtoTip, MAIN_DEPLOYMENT_ALTITUDE, drogue_drag_area, drogue_opening_time, main_drag_area, main_opening_time, main_altitude_time)  # N*m, moment due to parachute drag

    return M_body + M_rot + M_chute

# Calculate angular acceleration of rocket
def calc_alpha(t, alt, vx, vy, phi, omega, CGpos, CD_cylinder, rocket_radius, rocket_outer_diameter, total_length, total_mass, CGtoTip, MAIN_DEPLOYMENT_ALTITUDE, drogue_drag_area, drogue_opening_time, main_drag_area, main_opening_time, inertia, main_altitude_time):
    M_sum = calculate_moments(t, alt, omega, vx, vy, phi, CGpos, CD_cylinder, rocket_radius, rocket_outer_diameter, total_length, total_mass, CGtoTip, MAIN_DEPLOYMENT_ALTITUDE, drogue_drag_area, drogue_opening_time, main_drag_area, main_opening_time, main_altitude_time)

    alpha = M_sum / inertia  # rad/s^2, angular acceleration
    return alpha

def force_system(time, state, total_mass, CGpos, CD_cylinder, rocket_radius, rocket_outer_diameter, total_length, CGtoTip, MAIN_DEPLOYMENT_ALTITUDE, drogue_drag_area, drogue_opening_time, main_drag_area, main_opening_time, inertia, main_altitude_time):
    # unpack state vector
    phi, omega, velocity_x, velocity_y, altitude = state

    # Calculate forces
    Fg = Calculate_Force_Gravity(total_mass)
    Fs, main_altitude_time = Calculate_Force_Parachute_Drag(time, altitude, phi, omega, velocity_x, velocity_y, CGtoTip, MAIN_DEPLOYMENT_ALTITUDE, drogue_drag_area, drogue_opening_time, main_drag_area, main_opening_time, main_altitude_time)  # N, force due to parachute
    Fbd = Calculate_Force_Body_Drag(altitude, velocity_x, velocity_y, phi, total_length, total_mass, rocket_outer_diameter, CD_cylinder)

    F_total = Fg + Fs + Fbd

    acceleration = F_total / total_mass  # m/s^2
    ax, ay = acceleration

    # Calculate torque
    alpha = calc_alpha(time, altitude, velocity_x, velocity_y, phi, omega, CGpos, CD_cylinder, rocket_radius, rocket_outer_diameter, total_length, total_mass, CGtoTip, MAIN_DEPLOYMENT_ALTITUDE, drogue_drag_area, drogue_opening_time, main_drag_area, main_opening_time, inertia, main_altitude_time)  # rad/s^2

    return [omega, alpha, ax, ay, velocity_y]

def ground_check(time, state):
    return state[4]  # Event triggers when altitude reaches 0

def get_axial_load(phi, omega, mass_model, CGpos, total_mass, total_length, slice_pos, F_chute, F_body_drag, F_total):

    F_gs_axial =  [transform_gTor([0, Fg], phi)[0] for Fg in mass_model * -9.81] # N, gravitational force on each segment

    accel_axial = F_total[0]/total_mass  # m/s^2, axial acceleration of rocket
    F_chute_axial = F_chute[0] # N, axial component of parachute drag force
    Fbd_axial_dist = F_body_drag[0]/total_length  # N/m, distributed body drag force (q0)
    axial_inertial = accel_axial * mass_model
    axial = axial_inertial - F_gs_axial - Fbd_axial_dist * slice_pos
    axial[-1] -= F_chute_axial  # N, add parachute axial at top

    centripetal_force = np.abs(mass_model * pow(omega, 2) * (slice_pos - CGpos))
    axial += centripetal_force

    # axial_forces = np.abs(axial_total_forces) - np.abs(mass_model * axial_acceleration)
    # axial_forces = axial_total_forces - mass_model * axial_acceleration # + centripetal_force  # N, axial forces across rocket body
    return axial

def get_shear_bending(phi, alpha, total_length, CGpos, slice_pos, mass_model, mass_model_slices, total_mass, dy, F_parachute_drag, F_body_drag, F_total, M_rot):
    accel_shear = F_total[1]/total_mass  # m/s^2, shear acceleration of rocket
    F_gs_shear =  [transform_gTor([0, Fg], phi)[1] for Fg in mass_model * -9.81] # N, gravitational force on each segment
    Fbd_dist = F_body_drag[1]/total_length  # N/m, distributed body drag force (q0)

    shear_forces = np.zeros(len(slice_pos))
    for i, xk in enumerate(slice_pos):
        mass = mass_model[i]
        rot_mass = -np.sum(mass_model_slices[:i+1] * (CGpos - slice_pos[:i+1]))
        shear_forces[i] = accel_shear*mass + alpha*rot_mass
        shear_forces[i] -= F_gs_shear[i] + Fbd_dist * xk

    shear_forces[-1] -= F_parachute_drag[1]  # N, add parachute shear at top
    M_rot_dist = (M_rot/total_length) * slice_pos # Distribute rotational drag moment linearly along length

    bending_moments = np.cumsum(shear_forces) * dy + M_rot_dist

    return shear_forces, bending_moments

def plot(subplot, x_vectors, y_vectors, title, x_label, y_label, legends=['']):
    BACKGROUND_COLOR = "#e8e8e8"
    LINE_COLORS = ["#daaa00", "#8E6F3E", "#DDB945", "#CFB991"]

    subplot.set_facecolor(BACKGROUND_COLOR)
    if len(y_vectors) == len(x_vectors): # Single line
        y_vectors = [y_vectors]
    for subplot_count, y_vector in enumerate(y_vectors):
        subplot.plot(x_vectors, y_vector, color=LINE_COLORS[subplot_count], label=legends[subplot_count], linewidth=1)
    subplot.set_title(title)
    subplot.set_xlabel(x_label)
    subplot.set_ylabel(y_label)
    # subplot.grid()

    if legends != ['']:
        subplot.legend()

def simulate_recovery(show_plots = False):       
    
    # Plot styling
    plt.rcParams['lines.linewidth'] = 0.5
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Libre Franklin', 'Arial']
    plt.rcParams['axes3d.mouserotationstyle'] = 'azel'
    plt.rcParams["axes.grid"] = True
    plt.rcParams['legend.fontsize'] = 'small'

    PRIMARY_LINE_COLOR = "#DDB945"
    SECONDARY_LINE_COLOR = "#8E6F3E"

    
    loads.main()
    
    os.chdir(os.path.dirname(__file__))
    # print(f"os.getcwd: {os.getcwd()}")
        
    SIMULATION_DURATION = 800 # [seconds]
    time_step_size = 0.025 # [seconds]
    time_array = np.arange(0, SIMULATION_DURATION, time_step_size)


    if parameters.six_DoF_off_the_rail_velocity is None:
        apogee_altitude = parameters.one_DoF_estimated_apogee # [m] Apogee of rocket
        velocity_x_initial = 12 # [m/s]
    else:
        apogee_altitude = parameters.six_DoF_estimated_apogee # [m] Apogee of rocket
        velocity_x_initial = parameters.six_DoF_apogee_horizontal_velocity + 100000 # [m/s]    
    # apogee_altitude = parameters.six_DoF_estimated_apogee # [m] Apogee of rocket
    # velocity_x_initial = parameters.six_DoF_apogee_horizontal_velocity # [m/s]    

    dry_linear_density_array, length_along_rocket_linspace = get_mass_model(rocket_dict_dry, vehicle_parameters.parachute_mass)
    total_dry_mass = np.sum(dry_linear_density_array * (length_along_rocket_linspace[1] - length_along_rocket_linspace[0])) # [kg]

    dy = length_along_rocket_linspace[1] - length_along_rocket_linspace[0] # m, length of each slice in mass model
    dry_mass_model_slices = dry_linear_density_array * dy # [kg], mass of each slice in mass model
    dry_mass_model = np.cumsum(dry_mass_model_slices)  # [kg*m^2], mass model
    total_dry_mass = sum(dry_mass_model_slices)  # [kg]
    
    rocket_outer_diameter = parameters.tank_outer_diameter
    rocket_radius = rocket_outer_diameter/2  # [m]

    phi_init = 0 # [degrees] Initial orientation of rocket
    total_length = length_along_rocket_linspace[-1] # m
    length_vec = np.linspace(0, total_length, len(dry_mass_model_slices)) # [m] End position of each slice in mass model along rocket length, measured from aft
    slice_pos = length_vec - dy/2 # [m] Middle position of each slice in mass model along rocket length, measured from aft
    CGpos = np.sum(dry_mass_model_slices * slice_pos) / np.sum(dry_mass_model_slices) # CG position from aft (m), adjusted for no nosecone
    CGtoTip = total_length - CGpos  # CG position from tip (m)
    CGtoCenter = total_length/2 - CGpos  # CG position from center (m)

    # Moment of inertia calculation, parallel axis theorem
    I_slices = (1/12) * dry_mass_model_slices * (dy**2) # Moment of inertia of each slice about its own center (kg*m^2)
    I_parallel = dry_mass_model_slices * pow(slice_pos - CGpos, 2) # Mass * distance^2 term i n parallel axis theorem (kg*m^2)
    moment_of_inertia = np.sum(I_slices + I_parallel) + (1/4)*total_dry_mass*(rocket_radius**2) # kg*m^2, moment of inertia about CG

    CD_drogue = 2.2 # drag coefficient of drogue parachute
    CD_main = 0 # 2.2 # drag coefficient of main parachute

    CD_cylinder = 1.1

    initial_conditions = [phi_init, 0, velocity_x_initial, 0, apogee_altitude]  # [phi, omega, vx, vy, altitude]

    # Create chute area piecewise function
    drogue_radius = (168/2) * c.IN2M # rocket man 168 inch (bzb) parachute
    drogue_area = np.pi * pow(drogue_radius, 2)  # m^2
    drogue_drag_area = drogue_area * CD_drogue
    inflation_time_n = .86 # sec CHECK
    drogue_opening_time = inflation_time_n * (drogue_radius*2) / velocity_x_initial # seconds

    main_radius = 0.0000000000000000000000000000000000000001 * c.FT2M # m
    main_area = np.pi * pow(main_radius, 2)  # m^2
    main_drag_area = main_area * CD_main
    main_opening_time = 0.86 # seconds

    main_altitude_time = None # set this variable to the time at which the main chute opens
    MAIN_DEPLOYMENT_ALTITUDE = 0.01 # ft

    
    # idfk
    def ground_check_wrapper(time, state_vector, *unused_extra_arguments):
        return ground_check(time, state_vector)

    ground_check_wrapper.terminal = True
    ground_check_wrapper.direction = -1

    # Use solve_ivp to solve ODE system
    time_span = (0, SIMULATION_DURATION)
    solution = solve_ivp(force_system, 
                         time_span, 
                         initial_conditions, 
                         t_eval=time_array, 
                         events=ground_check_wrapper, 
                         method='RK45',
                         args=(total_dry_mass, 
                               CGpos, 
                               CD_cylinder, 
                               rocket_radius, 
                               rocket_outer_diameter, 
                               total_length, 
                               CGtoTip, 
                               MAIN_DEPLOYMENT_ALTITUDE, 
                               drogue_drag_area, 
                               drogue_opening_time, 
                               main_drag_area, 
                               main_opening_time,
                               moment_of_inertia,
                               main_altitude_time,), 
                        )

    # Unpack solution
    phi_array = solution.y[0]
    omega_array = solution.y[1]
    velocity_x_array = solution.y[2]
    velocity_y_array = solution.y[3]
    altitude_array = solution.y[4]

    if len(time_array) > len(solution.t):
        time_array = time_array[:len(solution.t)]
    else:
        raise ValueError
        # pass

    # Calculate forces across t_vec
    force_body_drag = np.array([Calculate_Force_Body_Drag(altitude_array[i], velocity_x_array[i], velocity_y_array[i], phi_array[i], total_length, total_dry_mass, rocket_outer_diameter, CD_cylinder) for i in range(len(velocity_x_array))])  # lbf, body drag force
    force_body_drag_magnitude = np.array([magnitude(f) for f in force_body_drag])
    force_body_drag_rocket_frame = np.array([transform_gTor(force_body_drag[i], phi_array[i]) for i in range(len(force_body_drag))])
    
    force_gravity = np.array([Calculate_Force_Gravity(total_dry_mass) for _ in time_array])
    force_gravity_magnitude = np.array([magnitude(f) for f in force_gravity])
    
    force_parachute_drag = np.array([Calculate_Force_Parachute_Drag(time_array[i], altitude_array[i], phi_array[i], omega_array[i], velocity_x_array[i], velocity_y_array[i], CGtoTip, MAIN_DEPLOYMENT_ALTITUDE, drogue_drag_area, drogue_opening_time, main_drag_area, main_opening_time, main_altitude_time)[0] for i in range(len(velocity_x_array))])  # lbf, parachute drag force
    
    for i in range(len(time_array)):
        _, main_altitude_time = Calculate_Force_Parachute_Drag(time_array[i], altitude_array[i], phi_array[i], omega_array[i], velocity_x_array[i], velocity_y_array[i], CGtoTip, MAIN_DEPLOYMENT_ALTITUDE, drogue_drag_area, drogue_opening_time, main_drag_area, main_opening_time, main_altitude_time)
        
        # print(f"i: {i}, main_altitude_time: {main_altitude_time}")
    if main_altitude_time is None:
        main_altitude_time = time_array[-1]
    
    
    force_parachute_drag_magnitude = np.array([magnitude(force) for force in force_parachute_drag])
    force_parachute_drag_rocket_frame = np.array([transform_gTor(force_parachute_drag[i], phi_array[i]) for i in range(len(force_parachute_drag))])  # N, parachute drag force in rocket frame

    # Sum forces and convert to Rocket frame
    force_total = force_gravity + force_body_drag + force_parachute_drag # N, total force
    force_total_rocket_frame = np.array([transform_gTor(force_total[i], phi_array[i]) for i in range(len(force_total))])  # lbf, total force in rocket frame

    # Calculate moments
    moment_due_to_rotational_drag = np.array([calculate_rotational_drag(altitude_array[i], omega_array[i], total_length, CD_cylinder, rocket_radius) for i in range(len(velocity_x_array))])  # N*m, moment due to rotational drag
    alphas = np.array([calc_alpha(time_array[i], altitude_array[i], velocity_x_array[i], velocity_y_array[i], phi_array[i], omega_array[i], CGpos, CD_cylinder, rocket_radius, rocket_outer_diameter, total_length, total_dry_mass, CGtoTip, MAIN_DEPLOYMENT_ALTITUDE, drogue_drag_area, drogue_opening_time, main_drag_area, main_opening_time, moment_of_inertia, main_altitude_time) for i in range(len(velocity_x_array))])  # rad/s^2, angular acceleration

    import time
    start_time = time.time()

    if show_plots == True:
        # make two subplots
        shear_bending_fig, (axial3d, shear3d, bending3d) = plt.subplots(1, 3, subplot_kw={'projection': '3d'})
        shear_bending_fig.suptitle("Internal Loads over Rocket Length over Time", fontsize=16)

    # 3D arrays for shear and bending
    shear_forces_over_time = []
    bending_moments_over_time = []
    axial_forces_over_time = []

    shear_bending_start_time = 0
    shear_bending_end_time = 10
    downsample_factor = 5
    
    axial_bending_start_time = 0
    axial_bending_end_time = 10
    
    time_array_shear_bending = np.arange(shear_bending_start_time, shear_bending_end_time, time_step_size)  # seconds
    time_array_axial = np.arange(shear_bending_start_time, shear_bending_end_time, time_step_size)  # seconds
    for shear_bending_time in time_array_shear_bending: # iterate over first __ seconds of flight
        mi = int(shear_bending_time/time_step_size) # Index for shear & bending
        # axial_t = index + main_altitude_time # Offset axial time to start at main deployment
        # axial_mi = int(axial_t/time_step_size)
        
        shear_bending_index = mi
        axial_index = mi
        # shear_bending_index = int(index)
        # axial_index = int(index + (main_altitude_time/time_step_size))
        
        shear, bending = get_shear_bending(phi_array[shear_bending_index], alphas[shear_bending_index], total_length, CGpos, slice_pos, dry_mass_model, dry_mass_model_slices, total_dry_mass, dy, force_parachute_drag_rocket_frame[shear_bending_index], force_body_drag_rocket_frame[shear_bending_index], force_total_rocket_frame[shear_bending_index], moment_due_to_rotational_drag[shear_bending_index])
        shear_forces_over_time.append(shear)  # [N] shear forces across rocket body
        bending_moments_over_time.append(bending)

        axial = get_axial_load(phi_array[axial_index], omega_array[axial_index], dry_mass_model, CGpos, total_dry_mass, total_length, slice_pos, force_parachute_drag_rocket_frame[axial_index], force_body_drag_rocket_frame[axial_index], force_total_rocket_frame[axial_index])
        axial_forces_over_time.append(axial)
    # plot surface on 3d plot
    length_3d, time_3d = np.meshgrid(np.array(length_vec) * c.M2FT, np.array(time_array_shear_bending), indexing='ij')

    if show_plots == True:
        bending3d.plot_surface(time_3d, length_3d, np.array(bending_moments_over_time).T*c.NM2LBI, cmap='viridis', rstride=downsample_factor, cstride=downsample_factor, alpha=0.9, edgecolor='none')
        bending3d.set_xlabel("Time (s)")
        bending3d.set_ylabel("Length from Aft (ft)")
        bending3d.set_zlabel("Bending Moment (ft-lbs)")
        bending3d.yaxis.set_major_locator(MultipleLocator(3))
        bending3d.set_title("Bending Moment")
        bending3d.invert_xaxis()

        shear3d.plot_surface(time_3d, length_3d, np.array(shear_forces_over_time).T*c.N2LBF, cmap='viridis', rstride=downsample_factor, cstride=downsample_factor, alpha=0.9, edgecolor='none')
        shear3d.set_xlabel("Time (s)")
        shear3d.set_ylabel("Length from Aft (ft)")
        shear3d.set_zlabel("Shear Force (lbf)")
        shear3d.yaxis.set_major_locator(MultipleLocator(3))
        shear3d.set_title("Shear Force")
        shear3d.invert_xaxis()

        # axial_time3d = np.array(time_array_axial)
        axial_time3d = time_3d

        axial3d.plot_surface(axial_time3d, length_3d, np.array(axial_forces_over_time).T * c.N2LBF, cmap='plasma', rstride=downsample_factor, cstride=downsample_factor, alpha=0.9, edgecolor='none')
        axial3d.set_xlabel("Time (s)")
        axial3d.set_ylabel("Length from Aft (ft)")
        axial3d.set_zlabel("Axial Load (lbf)")
        axial3d.set_title("Axial Load")
        axial3d.yaxis.set_major_locator(MultipleLocator(3))
        axial3d.invert_xaxis()

    # Find max. axial load, shear force & bending moment at each rocket point
    drogue_worst_bending_index = np.argmax(np.abs(np.array(bending_moments_over_time)), axis=0)[0]
    drogue_worst_axial_index = np.argmax(np.abs(np.array(axial_forces_over_time)), axis=0)[0]
    print("Drogue Worst Bending Time (s):", drogue_worst_bending_index * time_step_size)
    print(f"Drogue Worst Axial Time (s): {drogue_worst_axial_index * time_step_size:.5f}")
    
    
    # Get axial and bending at worst drogue bending time
    drogue_worst_bending_bending = bending_moments_over_time[drogue_worst_bending_index]
    drogue_worst_bending_shear = shear_forces_over_time[drogue_worst_bending_index]
    drogue_worst_bending_axial = get_axial_load(phi_array[drogue_worst_bending_index], omega_array[drogue_worst_bending_index],
                                        dry_mass_model,
                                        CGpos, 
                                        total_dry_mass, 
                                        total_length, 
                                        slice_pos, 
                                        force_parachute_drag_rocket_frame[drogue_worst_bending_index], 
                                        force_body_drag_rocket_frame[drogue_worst_bending_index], 
                                        force_total_rocket_frame[drogue_worst_bending_index])


    # Get axial and bending at worst drogue axial time
    drogue_worst_axial_bending = bending_moments_over_time[drogue_worst_axial_index]
    drogue_worst_axial_shear = shear_forces_over_time[drogue_worst_axial_index]
    drogue_worst_axial_axial = get_axial_load(phi_array[drogue_worst_axial_index], omega_array[drogue_worst_axial_index],
                                        dry_mass_model,
                                        CGpos, 
                                        total_dry_mass, 
                                        total_length, 
                                        slice_pos, 
                                        force_parachute_drag_rocket_frame[drogue_worst_axial_index], 
                                        force_body_drag_rocket_frame[drogue_worst_axial_index], 
                                        force_total_rocket_frame[drogue_worst_axial_index])





    main_worst_axial_index = np.argmax(np.abs(np.array(axial_forces_over_time)), axis=0)[0]
    print("Main Worst Axial Time (s):", (main_worst_axial_index * time_step_size) + main_altitude_time)
    # Get axial and bending at worst main axial time
    main_axial = axial_forces_over_time[main_worst_axial_index]
    main_shear, main_bending = get_shear_bending(phi_array[main_worst_axial_index],
                                        alphas[main_worst_axial_index],
                                        total_length, 
                                        CGpos, 
                                        slice_pos, 
                                        dry_mass_model, 
                                        dry_mass_model_slices, 
                                        total_dry_mass, 
                                        dy, 
                                        force_parachute_drag_rocket_frame[main_worst_axial_index], 
                                        force_body_drag_rocket_frame[main_worst_axial_index], 
                                        force_total_rocket_frame[main_worst_axial_index], 
                                        moment_due_to_rotational_drag[main_worst_axial_index])

    # Find worst shock load
    max_chute_load_time = np.argmax(force_parachute_drag_magnitude)
    max_chute_axial_load = transform_gTor(force_parachute_drag[max_chute_load_time], phi_array[max_chute_load_time])[0]

    mi = int(1/time_step_size)
    test_shear = shear_forces_over_time[mi]
    test_bending = bending_moments_over_time[mi]

    end_time = time.time()
    print(f"3D Graph Computation time: {end_time - start_time:.2f} seconds")
    print("-----------------------------------")
    print(f"Max axial load at recovery: {np.max(np.abs(drogue_worst_bending_axial))*c.N2LBF:.6f} lbf")
    print(f"Max shear force at recovery: {np.max(np.abs(drogue_worst_bending_shear))*c.N2LBF:.2f} lbf")
    print(f"Max bending moment at recovery: {np.max(np.abs(drogue_worst_bending_bending))*c.NM2LBI:.2f} ft-lbf")

    # Save to matlab
    # Converting to matlab file
    use_worst_axial = True
    if use_worst_axial:
        matlab_dict = {"axial_array": drogue_worst_axial_axial, "shear_array": drogue_worst_axial_shear, "bending_array": drogue_worst_axial_bending, "length_along_rocket_linspace": length_along_rocket_linspace} # Dictionary to save as .mat file
    else:
        matlab_dict = {"axial_array": drogue_worst_bending_axial, "shear_array": drogue_worst_bending_shear, "bending_array": drogue_worst_bending_bending, "length_along_rocket_linspace": length_along_rocket_linspace} # Dictionary to save as .mat file
    
    repository_root_path, _ = vehicle_parameters_functions.Get_Repository_Root_Path()
    savemat(repository_root_path / "SFD" / "rfd_outputs_recovery.mat", matlab_dict) # Save as .mat file for MATLAB)
    
    # Save internal loads to excel, for both drogue and main deployments
    import SFD.RDOF_2.output_formatter

    # Drogue
    #output_formatter.outputFormatter(shear=drogue_shear, bending=drogue_bending, axial=drogue_axial, shock_load=None, mass=total_mass, length=L, rocket_dict_aroldo=rocket_dict, sheet_name="Drogue")

    # Main
    #output_formatter.outputFormatter(shear=main_shear, bending=main_bending, axial=main_axial, shock_load=max_chute_axial_load, mass=total_mass, length=L, rocket_dict_aroldo=rocket_dict, sheet_name="Main")

    if show_plots == True:
        plt.show() # Show 3D plots

        # Create figure with 3 subplots: Axial, shear & bending
        internal_fig, (drogue_shear_plot, drogue_bending_plot, drogue_axial_plot, main_shear_plot, main_bending_plot, main_axial_plot) = plt.subplots(1,6, figsize=(12, 4))
        internal_fig.suptitle("Worst-Case Internal Rocket Loads during Recovery", fontsize=16)
        WINDOW_COLOR = "white"
        internal_fig.patch.set_facecolor(WINDOW_COLOR)

        # plot(drogue_shear_plot, length_vec * c.M2FT, np.array(drogue_worst_bending_shear)*c.N2LBF, 'Drogue - Shear Force', 'Length from Aft (ft)', 'Shear Force (lbf)')
        # plot(drogue_bending_plot, length_vec * c.M2FT, np.array(drogue_worst_bending_bending)*c.NM2LBI, 'Drogue - Bending Moment', 'Length from Aft (ft)', 'Bending Moment (ft-lbs)')
        # plot(drogue_axial_plot, length_vec * c.M2FT, np.array(drogue_worst_bending_axial)*c.N2LBF, 'Drogue - Axial Load', 'Length from Aft (ft)', 'Axial Load (lbf)')
        plot(drogue_shear_plot, length_vec * c.M2FT, np.array(drogue_worst_axial_shear)*c.N2LBF, 'Drogue - Shear Force', 'Length from Aft (ft)', 'Shear Force (lbf)')
        plot(drogue_bending_plot, length_vec * c.M2FT, np.array(drogue_worst_axial_bending)*c.NM2LBI, 'Drogue - Bending Moment', 'Length from Aft (ft)', 'Bending Moment (ft-lbs)')
        plot(drogue_axial_plot, length_vec * c.M2FT, np.array(drogue_worst_axial_axial)*c.N2LBF, 'Drogue - Axial Load', 'Length from Aft (ft)', 'Axial Load (lbf)')
        
        plot(main_shear_plot, length_vec * c.M2FT, np.array(main_shear)*c.N2LBF, 'Main - Shear Force', 'Length from Aft (ft)', 'Shear Force (lbf)')
        plot(main_bending_plot, length_vec * c.M2FT, np.array(main_bending)*c.NM2LBI, 'Main - Bending Moment', 'Length from Aft (ft)', 'Bending Moment (ft-lbs)')
        plot(main_axial_plot, length_vec * c.M2FT, np.array(main_axial)*c.N2LBF, 'Main - Axial Load', 'Length from Aft (ft)', 'Axial Load (lbf)')
        
        plt.tight_layout()
        plt.show()

        altitude_array *= c.M2FT
        velocity_x_array *= c.M2FT
        velocity_y_array *= c.M2FT

        # Plot a 2x4 grid of subplots:
        fig, axs = plt.subplots(2, 4, figsize=(12, 8))
        fig.suptitle(f"RDOF - {drogue_radius*2*c.M2FT:.1f} ft. drogue, {main_radius*2*c.M2FT:.1f} ft. main", fontsize=16)
        fig.patch.set_facecolor(WINDOW_COLOR)

        # Subplot 1: Velocity x and y vs time
        plot(axs[0,0], time_array, [velocity_x_array, velocity_y_array], 'Velocity vs Time', 'Time (s)', 'Velocity (ft/s)', legends=['Vx', 'Vy'])

        # Subplot 2: Individual forces vs time
        plot(axs[0,1], time_array, [force_gravity_magnitude*c.N2LBF, force_body_drag_magnitude*c.N2LBF, force_parachute_drag_magnitude*c.N2LBF], 'Global Forces vs Time', 'Time (s)', 'Force (lbf)', legends=['Weight', 'Body Drag', 'Chute Drag'])

        # Subplot 3: Phi and velocity angle vs time
        phi_vertical = np.pi/2 - phi_array
        plot(axs[1,0], time_array, [phi_vertical * 180/
                            np.pi], 'Rocket Angle vs Time', 'Time (s)', 'Angle from Vertical (degrees)', legends=['Rocket Angle'])

        # Subplot 4: Axial Load across rocket body
        plot(axs[1,1], length_vec * c.M2FT, test_shear * c.N2LBF, 'Shear Force', 'Length from Aft (ft)', 'Shear Force (lbf)')
        # Subplot 5: Rocket frame forces vs time
        plot(axs[0,2], time_array, [force_total_rocket_frame[:, 0]*c.N2LBF, force_total_rocket_frame[:, 1]*c.N2LBF], 'Rocket Frame Forces vs Time', 'Time (s)', 'Force (lbf)', legends=['Axial', 'Shear'])
        # Subplot 6: Internal bending moments
        plot(axs[1,2], length_vec * c.M2FT, test_bending * c.NM2LBI, 'Bending Moment', 'Length from Aft (ft)', 'Bending Moment (ft-lbs)')
        # Subplot 7: Altitude vs time
        plot(axs[0,3], time_array, altitude_array, 'Altitude vs Time', 'Time (s)', 'Altitude (ft)')

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

            string = f"t = {time_array[i]:.1f}s\nAltitude = {altitude_array[i]:.1f}ft"
            title.set_text(string)

            rocket_angle = phi_vertical[i] # degrees
            rocket_angle_rad = rocket_angle

            # update line to point in rocket_angle_rad direction, and length of line is CGtoTip
            line.set_data([rocket_angle_rad, rocket_angle_rad], [0, CGtoTip])
            bottomLine.set_data([rocket_angle_rad+np.pi, rocket_angle_rad+np.pi], [0, CGpos])

            return line, bottomLine, title

        anim = FuncAnimation(fig, update, frames=len(time_array)//anim_speed_factor, blit=True, interval=1, repeat=True)
        
        # Save animation as MP4 file
        # anim.save('RDOF_animation.mp4', writer='ffmpeg', fps=30, dpi=300)

        # # Save the figure as a PNG file
        # fig.savefig('RDOF.png', dpi=300, bbox_inches='tight')

        # Find ground impact speed
        # plt.figure()
        # plt.plot(altitude_array)
        
        # Show animation
        # mng = plt.get_current_fig_manager()
        # mng.window.showMaximized()
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()
        
    groundIndex = np.where(altitude_array <= 5)
    if groundIndex[0].size > 0:
        ind = groundIndex[0][0]
        descent_time = ind*time_step_size
        print(f"descent time: {descent_time:.2f} seconds")
        velocity = magnitude([velocity_x_array[ind], velocity_y_array[ind]])
        print(f"Rocket velocity at ground: {velocity:.2f} ft/s")
    else:
        print("Rocket doesn't reach ground given timeframe/conditions.")


def main():
    simulate_recovery()

if __name__ == "__main__":
    rerun_everything = False
    
    # change this to somehow work to just up to this file
    if rerun_everything:
        vehicle_main.vehicle_analysis()
    else:
        pass
        # parameters = csv_to_dataclass(parameters_csv_filepath)    
        
    simulate_recovery(show_plots = True)