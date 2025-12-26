import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
import main as m
import vehicle_parameters as v

#ATMOSPHERE_DATA = pd.read_csv("atmosphere.csv")
#ATMOSPHERE_DATA = np.array(ATMOSPHERE_DATA)

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
    reference_area = (np.pi * (tankOD) ** 2 / 4) + finNumber * finHeight * c.FIN_THICKNESS # [m^2] reference area of the rocket
    mass = wetMass  # [kg] initial mass of the rocket

    cD = 0.5
    ascent_drag_coefficient = cD * (totalLength / 6.35) * (tankOD / 0.203)

    # Initial Conditions
    altitude = c.INDIANA_ALTITUDE  # [m] initial altitude of the rocket
    velocity = 0  # [m/s] initial velocity of the rocket
    acceleration = ((jetThrust + (exitPressure - ATMOSPHERE_DATA[(0, 1)]) * exitArea) - (c.GRAVITY * wetMass)) / wetMass  # [m/s] initial acceleration of the rocket
    time = 0  # [s] initial time of the rocket
    dt = 0.0001  # [s] time step of the rocket

    # Array Initialization:
    altitude_array = [altitude]
    velocity_array = [velocity]
    acceleration_array = [acceleration]
    time_array = [time]
    
    totalImpulse = 0  # Initialize total impulse

    count = 1
    while (velocity >= 0) or (acceleration > 0):

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
        
        if time < burnTime:
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
        da_TWR = thrust / weight  # acceleration equation of motion

        acceleration = (thrust - drag_force - weight) / mass  # acceleration equation of motion
        acceleration_array.append(acceleration)

        velocity += acceleration * dt  # velocity integration
        velocity_array.append(velocity)

        altitude += velocity * dt  # position integration
        altitude_array.append(altitude)

        time = time + dt  # time is inevitable
        time_array.append(time)
        
        count += 1
        
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

def mdot_based_on_time(mDotTotal):
    angle_history = m.angle_history
    time_history = m.time_history
    actuation_time = m.time
    mdot_history = []
    csv_path = os.path.join(os.path.dirname(__file__), "anglevscv.csv")
    angle_data = pd.read_csv(csv_path)
    angles = angle_data.iloc[:, 0].to_numpy()
    cv_frac = angle_data.iloc[:, 1].to_numpy()
    f_linear = interp1d(angles, cv_frac, kind='linear', bounds_error=False, fill_value=(0.0, 1.0))

    if time_history[2] - time_history[1] != 0.0001:
        print("Ensure time step of 0.0001 in actuator sizing code!")
        return None
    for time, angle in zip(time_history, angle_history):
            cv_fraction = float(f_linear(angle))
            mdot = mDotTotal * cv_fraction
            mdot_history.append(mdot)
            
    mdot_history = np.array(mdot_history)
    plt.plot(time_history, mdot_history * c.KG2LBM)
    plt.xlabel('Time [s]')
    plt.ylabel('mdot [lbm/s]')
    plt.show()
    plt.plot(angle_history, mdot_history * c.KG2LBM)
    plt.xlabel('Valve Angle [degrees]')
    plt.ylabel('mdot [lbm/s]')
    plt.show()
    return None

mDotTotal = v.parameters.total_mass_flow_rate
mdot_based_on_time(mDotTotal)