import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import main as m
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
import vehicle_parameters as v

csv_path_atmosphere = os.path.join(os.path.dirname(__file__), "atmosphere.csv")
ATMOSPHERE_DATA = pd.read_csv(csv_path_atmosphere)
ATMOSPHERE_DATA = np.array(ATMOSPHERE_DATA)

def calculate_trajectory(
    wetMass,
    mDotTotal,
    jetThrust,
    tankOD,
    finNumber,
    finHeight,
    exitArea,
    exitPressure,
    burnTime,
    totalLength,
    plots,
):
    """
    _summary_

    Parameters
    ----------
    wetMass : float
        Wet mass of the rocket [kg].
    mDotTotal : float
        Maximum total mass flow rate of the engine [kg/s].
    jetThrust : float
        Engine thrust [N].
    tankOD : float
        Outer diameter of the tank [m].
    finNumber : int
        Number of fins [-].
    finHeight : float
        Fin semi-span [m].
    exitArea : float
        Exit area of the nozzle [m^2].
    exitPressure : float
        Exit pressure of the nozzle [Pa].
    burnTime : float
        Burn time of the engine [s].
    totalLength : float
        Total Length of Rocket [m].
    plots : bool
        Boolean for plotting, 1 = on, 0 = off [-].

    Returns
    -------
    altitude : float
        Final altitude of the rocket [m].
    maxMach : float
        Maximum Mach number of the rocket [-].
    maxAccel : float
        Maximum acceleration of the rocket [m/s^2].
    exitVelo : float
        Exit velocity of the rocket [m/s].
    totalImpulse : float
        Total impulse of the rocket [Ns].
    """

    # Rocket Properties
    reference_area = (np.pi * (tankOD) ** 2 / 4) + finNumber * finHeight * 1/2 * c.IN2M # [m^2] reference area of the rocket, 1/2 in fin thickness
    mass = wetMass  # [kg] initial mass of the rocket

    cD = 0.5
    ascent_drag_coefficient = cD * (totalLength / 6.35) * (tankOD / 0.203)

    # Initial Conditions
    altitude = c.INDIANA_ALTITUDE  # [m] initial altitude of the rocket
    velocity = 0  # [m/s] initial velocity of the rocket
    acceleration = 0  # [m/s] initial acceleration of the rocket
    time = 0  # [s] initial time of the rocket
    dt = 0.0001  # [s] time step of the rocket

    # Array Initialization:
    altitude_array = [altitude]
    velocity_array = [velocity]
    acceleration_array = [acceleration]
    time_array = [time]
    
    totalImpulse = 0  # Initialize total impulse

    #mdot_history, main_thrust_history, time_history = mdot_based_on_time(mDotTotal, jetThrust)
    actuation_time = 0.5 # hardcoded actuation time
    mdot_history, main_thrust_history, time_history = mdot_hardcoded_actuation(mDotTotal, actuation_time, jetThrust)

    count = 1
    while (velocity >= 0) or (acceleration >= 0):

        altitude_index = int(altitude // 10)  # Divide altitude by 10 to find index

        if altitude_index < 0:
            pressure = ATMOSPHERE_DATA[(0, 1)]
            rho = ATMOSPHERE_DATA[(0, 2)]  # Return first row if below range
        elif altitude_index >= len(ATMOSPHERE_DATA):
            pressure = ATMOSPHERE_DATA[(-1, 1)]
            rho = ATMOSPHERE_DATA[(-1, 2)]  # Return last row if above range
        else:
            pressure = ATMOSPHERE_DATA[(altitude_index, 1)]
            rho = ATMOSPHERE_DATA[(altitude_index, 2)]
        
        if time <= time_history[-1] and time != 0:
            mass = mass - mdot_history[int(time / dt)] * dt  # [kg] mass of the rocket
            thrust = (
                main_thrust_history[int(time / dt)] + (exitPressure - pressure) * exitArea
            )  # [N] force of thrust, accounting for pressure thrust

            totalImpulse += thrust * dt  # Accumulate impulse
            if mass < 0:
                raise ValueError("youre a dumbass, CHECK IF TANK DATA IS FAKE OR NOT!!!!!!!!!!!!!")
        elif time < burnTime and time > time_history[-1]:
            mass = mass - mDotTotal * dt  # [kg] mass of the rocket
            thrust = (
                jetThrust + (exitPressure - pressure) * exitArea
            )  # [N] force of thrust, accounting for pressure thrust

            totalImpulse += thrust * dt  # Accumulate impulse
            if mass < 0:
                raise ValueError("youre a dumbass, CHECK IF TANK DATA IS FAKE OR NOT!!!!!!!!!!!!!")
        else:
            thrust = 0  # [N] total thrust of the rocket

        drag_force = (
            0.5 * rho * velocity**2 * ascent_drag_coefficient * reference_area
        )  # [N] force of drag
        weight = c.GRAVITY * mass  # [N] downward force due to gravity

        acceleration = (thrust - drag_force - weight) / mass  # acceleration equation of motion
        acceleration_array.append(acceleration)

        velocity += acceleration * dt  # velocity integration
        altitude += velocity * dt  # position integration
        time = time + dt  # time is inevitable

        if velocity < 0 and altitude <= altitude_array[0]:
            velocity = 0
            altitude = altitude_array[0]
        velocity_array.append(velocity)
        altitude_array.append(altitude)
        time_array.append(time)
        count += 1
    print(acceleration_array[:50])
    # Find the closest altitude to the TRAILER_RAIL_HEIGHT

    off_the_rail_velocity, off_the_rail_acceleration, off_the_rail_time = 0, 0, 0
    for time_step_index in range(len(altitude_array)):
        if altitude_array[time_step_index] >= c.TRAILER_RAIL_HEIGHT + c.INDIANA_ALTITUDE:
            off_the_rail_velocity = velocity_array[time_step_index]
            off_the_rail_acceleration = acceleration_array[time_step_index]
            off_the_rail_time = (time_step_index/len(altitude_array)) * time_array[-1]
            break

    estimated_apogee = altitude * 0.651 # what the fuck is this for

    if plots == True:
        plt.plot(
                 time_array, 
                 np.array(acceleration_array, dtype=float) * c.M2FT,
                 label="Acceleration (any direction) [ft/s^2]",
                 )
      
        plt.plot(
                 time_array, 
                 np.array(velocity_array, dtype=float) * c.M2FT,
                 label="Velocity (any direction) [ft/s]",
                 )
   
        
        plt.plot(time_array,
                 np.array(altitude_array, dtype=float) * c.M2FT * 0.651,
                label="Height [ft]"
                )
        plt.title("Height, Velocity, Acceleration v. Time")
        plt.legend()
        
        plt.ylabel("Height [ft], Velocity (any direction) [ft/s], Acceleration (any direction) [ft/s^2]")
        plt.xlabel("Time [s]")
        plt.grid()
        plt.show()


    max_velocity = velocity_array[np.argmax(np.abs(velocity_array))] # account for positive or negative max
    max_acceleration = acceleration_array[np.argmax(np.abs(acceleration_array))] # account for positive or negative max
    
    
    return (
        float(estimated_apogee),
        max_acceleration,
        max_velocity,
        float(off_the_rail_velocity),
        float(off_the_rail_acceleration),
        float(totalImpulse),
        float(off_the_rail_time),
    )

def calculate_trajectory_variable_burn(
    wetMass,
    mDotTotal,
    jetThrust,
    tankOD,
    finNumber,
    finHeight,
    exitArea,
    exitPressure,
    totalLength,
    propellant_mass,
    plots,
):
    """
    _summary_

    Parameters
    ----------
    wetMass : float
        Wet mass of the rocket [kg].
    mDotTotal : float
        Maximum total mass flow rate of the engine [kg/s].
    jetThrust : float
        Engine thrust [N].
    tankOD : float
        Outer diameter of the tank [m].
    finNumber : int
        Number of fins [-].
    finHeight : float
        Fin semi-span [m].
    exitArea : float
        Exit area of the nozzle [m^2].
    exitPressure : float
        Exit pressure of the nozzle [Pa].
    burnTime : float
        Burn time of the engine [s].
    totalLength : float
        Total Length of Rocket [m].
    plots : bool
        Boolean for plotting, 1 = on, 0 = off [-].

    Returns
    -------
    altitude : float
        Final altitude of the rocket [m].
    maxMach : float
        Maximum Mach number of the rocket [-].
    maxAccel : float
        Maximum acceleration of the rocket [m/s^2].
    exitVelo : float
        Exit velocity of the rocket [m/s].
    totalImpulse : float
        Total impulse of the rocket [Ns].
    """

    # Rocket Properties
    reference_area = (np.pi * (tankOD) ** 2 / 4) + finNumber * finHeight * 1/2 * c.IN2M # [m^2] reference area of the rocket, 1/2 in fin thickness
    mass = wetMass  # [kg] initial mass of the rocket

    cD = 0.5
    ascent_drag_coefficient = cD * (totalLength / 6.35) * (tankOD / 0.203)

    # Initial Conditions
    altitude = c.INDIANA_ALTITUDE  # [m] initial altitude of the rocket
    velocity = 0  # [m/s] initial velocity of the rocket
    acceleration = 0  # [m/s] initial acceleration of the rocket
    time = 0  # [s] initial time of the rocket
    dt = 0.0001  # [s] time step of the rocket

    # Array Initialization:
    altitude_array = [altitude]
    velocity_array = [velocity]
    acceleration_array = [acceleration]
    time_array = [time]
    
    totalImpulse = 0  # Initialize total impulse

    #mdot_history, main_thrust_history, time_history = mdot_based_on_time(mDotTotal, jetThrust)
    actuation_time = 0.5 # hardcoded actuation time
    mdot_history, main_thrust_history, time_history = mdot_hardcoded_actuation(mDotTotal, actuation_time, jetThrust)
    burnTime = 0 #### Calculate the burn time
    
    for mdot in mdot_history:
            mass_lost = mdot * dt
            burnTime += dt
            propellant_mass -= mass_lost
    while propellant_mass > 0:
        mass_lost = mdot_history[-1] * dt
        burnTime += dt
        propellant_mass -= mass_lost

    count = 1
    while (velocity >= 0) or (acceleration >= 0):

        altitude_index = int(altitude // 10)  # Divide altitude by 10 to find index

        if altitude_index < 0:
            pressure = ATMOSPHERE_DATA[(0, 1)]
            rho = ATMOSPHERE_DATA[(0, 2)]  # Return first row if below range
        elif altitude_index >= len(ATMOSPHERE_DATA):
            pressure = ATMOSPHERE_DATA[(-1, 1)]
            rho = ATMOSPHERE_DATA[(-1, 2)]  # Return last row if above range
        else:
            pressure = ATMOSPHERE_DATA[(altitude_index, 1)]
            rho = ATMOSPHERE_DATA[(altitude_index, 2)]
        
        if time <= time_history[-1] and time != 0:
            mass = mass - mdot_history[int(time / dt)] * dt  # [kg] mass of the rocket
            thrust = (
                main_thrust_history[int(time / dt)] + (exitPressure - pressure) * exitArea
            )  # [N] force of thrust, accounting for pressure thrust

            totalImpulse += thrust * dt  # Accumulate impulse
            if mass < 0:
                raise ValueError("youre a dumbass, CHECK IF TANK DATA IS FAKE OR NOT!!!!!!!!!!!!!")
        elif time < burnTime and time > time_history[-1]:
            mass = mass - mDotTotal * dt  # [kg] mass of the rocket
            thrust = (
                jetThrust + (exitPressure - pressure) * exitArea
            )  # [N] force of thrust, accounting for pressure thrust

            totalImpulse += thrust * dt  # Accumulate impulse
            if mass < 0:
                raise ValueError("youre a dumbass, CHECK IF TANK DATA IS FAKE OR NOT!!!!!!!!!!!!!")
        else:
            thrust = 0  # [N] total thrust of the rocket

        drag_force = (
            0.5 * rho * velocity**2 * ascent_drag_coefficient * reference_area
        )  # [N] force of drag
        weight = c.GRAVITY * mass  # [N] downward force due to gravity

        acceleration = (thrust - drag_force - weight) / mass  # acceleration equation of motion
        acceleration_array.append(acceleration)

        velocity += acceleration * dt  # velocity integration
        altitude += velocity * dt  # position integration
        time = time + dt  # time is inevitable

        if velocity < 0 and altitude <= altitude_array[0]:
            velocity = 0
            altitude = altitude_array[0]
        velocity_array.append(velocity)
        altitude_array.append(altitude)
        time_array.append(time)
        count += 1
    print(acceleration_array[:50])
    # Find the closest altitude to the TRAILER_RAIL_HEIGHT

    off_the_rail_velocity, off_the_rail_acceleration, off_the_rail_time = 0, 0, 0
    for time_step_index in range(len(altitude_array)):
        if altitude_array[time_step_index] >= c.TRAILER_RAIL_HEIGHT + c.INDIANA_ALTITUDE:
            off_the_rail_velocity = velocity_array[time_step_index]
            off_the_rail_acceleration = acceleration_array[time_step_index]
            off_the_rail_time = (time_step_index/len(altitude_array)) * time_array[-1]
            break

    estimated_apogee = altitude * 0.651 # what the fuck is this for

    if plots == True:
        plt.plot(
                 time_array, 
                 np.array(acceleration_array, dtype=float) * c.M2FT,
                 label="Acceleration (any direction) [ft/s^2]",
                 )
      
        plt.plot(
                 time_array, 
                 np.array(velocity_array, dtype=float) * c.M2FT,
                 label="Velocity (any direction) [ft/s]",
                 )
   
        
        plt.plot(time_array,
                 np.array(altitude_array, dtype=float) * c.M2FT * 0.651,
                label="Height [ft]"
                )
        plt.title("Height, Velocity, Acceleration v. Time")
        plt.legend()
        
        plt.ylabel("Height [ft], Velocity (any direction) [ft/s], Acceleration (any direction) [ft/s^2]")
        plt.xlabel("Time [s]")
        plt.grid()
        plt.show()


    max_velocity = velocity_array[np.argmax(np.abs(velocity_array))] # account for positive or negative max
    max_acceleration = acceleration_array[np.argmax(np.abs(acceleration_array))] # account for positive or negative max
    
    
    return (
        float(estimated_apogee),
        max_acceleration,
        max_velocity,
        float(off_the_rail_velocity),
        float(off_the_rail_acceleration),
        float(totalImpulse),
        float(off_the_rail_time),
    )

def mdot_based_on_time(mDotTotal, jetThrust):
    dt = 0.0001
    angle_history = m.angle_history
    time_history = m.time_history
    main_thrust_history = []
    mdot_history = []
    csv_path = os.path.join(os.path.dirname(__file__), "anglevscv.csv")
    angle_data = pd.read_csv(csv_path)
    angles = angle_data.iloc[:, 0].to_numpy()
    cv_frac = angle_data.iloc[:, 1].to_numpy()
    f_linear = interp1d(angles, cv_frac, kind='linear', bounds_error=False, fill_value=(0.0, 1.0))

    if time_history[2] - time_history[1] != dt:
        print("Ensure time step of 0.0001 in actuator sizing code!")
        return None
    for angle in angle_history:
            cv_fraction = float(f_linear(angle))
            mdot = mDotTotal * cv_fraction
            mdot_history.append(mdot)
            main_thrust = (mdot / mDotTotal) * jetThrust
            main_thrust_history.append(main_thrust)
    main_thrust_history = np.array(main_thrust_history)
    mdot_history = np.array(mdot_history)
    plt.plot(time_history, mdot_history * c.KG2LBM)
    plt.xlabel('Time [s]')
    plt.ylabel('mdot [lbm/s]')
    plt.show()
    plt.plot(angle_history, mdot_history * c.KG2LBM)
    plt.xlabel('Valve Angle [degrees]')
    plt.ylabel('mdot [lbm/s]')
    plt.show()
    plt.plot(time_history, main_thrust_history)
    plt.xlabel('Time [s]')
    plt.ylabel('Main Thrust [N]')
    plt.show()
    print(f"Main Thrust at kast: {main_thrust}")
    return mdot_history, main_thrust_history, time_history

def mdot_hardcoded_actuation(mDotTotal, actuation_time, jetThrust):
    dt = 0.0001
    time_history = np.arange(0, actuation_time + dt, dt)
    mdot_history = mDotTotal * (time_history / actuation_time)
    main_thrust_history = (mdot_history / mDotTotal) * jetThrust
    return mdot_history, main_thrust_history, time_history

def actuation_time_vs_rail_exit():
    max_actuation_time = 1 # [second]
    for actuation_time in np.arange(0, max_actuation_time, 0.0001):
        pass

        
#############################
#######  PARAMETERS  ########
#############################
wetMass = v.vehicle_wet_mass
mDotTotal = v.parameters.total_mass_flow_rate
jetThrust = v.parameters.jet_thrust
tankOD = v.parameters.tank_outer_diameter
finNumber =  v.number_of_fins
finHeight = 12 * c.IN2M
exitArea = np.pi * 4**2 / 4 * c.IN22M2
exitPressure = v.parameters.exit_pressure
burnTime = v.parameters.burn_time
totalLength = v.rocket_length
propellant_mass = v.parameters.total_propellant_mass
plots = 1

estimated_apogee, max_acceleration, max_velocity, off_the_rail_velocity, off_the_rail_acceleration, totalImpulse, off_the_rail_time = calculate_trajectory_variable_burn(
    wetMass,
    mDotTotal,
    jetThrust,
    tankOD,
    finNumber,
    finHeight,
    exitArea,
    exitPressure,
    totalLength,
    propellant_mass,
    plots)

print(estimated_apogee, max_acceleration, max_velocity, off_the_rail_velocity, off_the_rail_acceleration, totalImpulse, off_the_rail_time)