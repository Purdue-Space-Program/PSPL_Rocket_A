import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import constants as c
import vehicle_parameters
parameters, wet_mass_distribution, dry_mass_distribution = vehicle_parameters.main()
import vehicle_parameters_functions
import vehicle_main
import parseWind
import sfd



# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------
def max_q_or_off_the_rail(plot_on = False, print_inputs = False):
    # Select location to analyze
    location = "max_q" # Change to "max_q" or "off_the_rail"
    # ------------------------------------------------------------------------------

    # Inputs
    diameter = parameters.tube_outer_diameter # [m]
    thrust = parameters.jet_thrust # [N]
    total_length = parameters.total_length # [m]

    air_density = 1.225 # [kg / m^3]
    max_q_wind_gust = parseWind.percentile_75_wind_gust_speed # [m / s]
    off_the_rail_rail_whip = 5 # [m / s] about 11 mph

    if location == "max_q":
        wind_gust = max_q_wind_gust # [m/s] Wind gust at max q
        CoM = parameters.dry_COM_location_from_bottom # [m] Center of gravity from bottom of rocket
        linear_density_array, length_along_rocket_linspace = sfd.mass_model(vehicle_parameters.rocket_dict_dry)
        total_mass = parameters.dry_mass # [kg] Mass at max q
        if parameters.six_DoF_max_Q_velocity is None:
            velocity = parameters.one_DoF_max_Q_velocity # [m / s] Velocity at max q
            acceleration = parameters.one_DoF_max_Q_acceleration # [m / s^2] Acceleration at max q
            mach = parameters.one_DoF_max_Q_mach # [dimensionless] Mach number at max q
        else:
            velocity = parameters.six_DoF_max_Q_velocity # [m / s] Velocity at max q
            acceleration = parameters.six_DoF_max_Q_acceleration # [m / s^2] Acceleration at max q
            mach = parameters.six_DoF_max_Q_mach # [dimensionless] Mach number at max q
    
    elif location == "off_the_rail":
        wind_gust = off_the_rail_rail_whip # [m/s] Rail whip off the rail
        CoM = parameters.wet_COM_location_from_bottom # [m] Center of gravity from bottom of rocket
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        ######################################################################################################################################################
        ######################################################################################################################################################
        ######################################################################################################################################################
        ######################################################################################################################################################
        ######################################################################################################################################################
        ############################################# Replace Rocket_Dict_Wet And Dry With Wet_Mass Distribution!!!!!! #######################################
        ######################################################################################################################################################
        ######################################################################################################################################################
        ######################################################################################################################################################
        ######################################################################################################################################################
        ######################################################################################################################################################
        ######################################################################################################################################################
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        linear_density_array, length_along_rocket_linspace = sfd.mass_model(vehicle_parameters.rocket_dict_wet)
        total_mass = parameters.wet_mass # [kg] Mass off the rail
        if parameters.six_DoF_off_the_rail_velocity is None:
            velocity = parameters.one_DoF_off_the_rail_velocity # [m / s] Velocity off the rail
            acceleration = parameters.one_DoF_off_the_rail_acceleration # [m / s^2] Acceleration off the rail
            mach = velocity / c.SPEED_OF_SOUND # [Mach] Speed of sound near sea level
        else:
            velocity = parameters.six_DoF_off_the_rail_velocity # [m / s] Velocity off the rail
            acceleration = parameters.six_DoF_off_the_rail_acceleration # [m / s^2] Acceleration off the rail
            mach = velocity / c.SPEED_OF_SOUND # [Mach] Speed of sound near sea level
    else:
        raise ValueError("Invalid location selected. Choose 'max_q' or 'off_the_rail'.")

    total_mass = np.sum(linear_density_array * (length_along_rocket_linspace[1] - length_along_rocket_linspace[0])) # [kg]

    dx = length_along_rocket_linspace[1] - length_along_rocket_linspace[0]  # [m]

    # Fins
    root_chord = parameters.root_chord # [m]
    tip_chord = parameters.tip_chord # [m]
    sweep_length = parameters.sweep_length # [m]
    fin_height = 7.25 * c.IN2M # [m] what is this, why is it not just root_chord
    numFins = 3 # [m]
    fin_top = parameters.fin_top # [m]
    noseconeToFin = total_length - fin_top # [m]

    # Calculated inputs
    Q = sfd.calcQ(air_density, velocity) # [N/m^2] Dynamic pressure
    AOA = sfd.calcAOA(wind_gust, velocity)
    S = sfd.calcS(diameter) # [m^2] Cross sectional area
    # ------------------------------------------------------------------------------

    # Calculated values
    finCP = sfd.calcFinCP(root_chord, tip_chord, sweep_length, fin_height, total_length, noseconeToFin) # Fin center of pressure
    noseCP = sfd.calcNoseCP(wet_mass_distribution.nosecone.length, total_length) # Nose center of pressure
    finSD = sfd.calcFinSD(root_chord, tip_chord, sweep_length, fin_height, numFins, diameter) # Fin stability derivative
    machCoeff = sfd.calcMachCoeff(1, mach) # Mach coefficient
    noseSD = sfd.calcNoseSD(CoM, noseCP, diameter) # Nose stability derivative
    noseLift = sfd.calcLift(Q, S, AOA, noseSD) # Nose lift
    finLift = sfd.calcLift(Q, S, AOA, finSD) # Fin lift
    inertia = sfd.calcRotationalInertia(linear_density_array, length_along_rocket_linspace, CoM) # Rotational inertia
    ay = sfd.calcLateralAcceleration(noseLift, finLift, total_mass) # Lateral acceleration
    r = sfd.calcAngularAcceleration(noseLift, finLift, noseCP, finCP, inertia, CoM) # Angular acceleration
    shear_array = np.array(sfd.calcShear(noseLift, finLift, noseCP, finCP, ay, linear_density_array, length_along_rocket_linspace, r, CoM)) # Shear force array
    bending_array = np.array(sfd.calcBending(shear_array, length_along_rocket_linspace)) # Bending moment array
    axial_array = np.array(sfd.calcAxial(thrust, acceleration, linear_density_array, length_along_rocket_linspace, air_density, 0.65, S, velocity)) # Axial forces array, For medium size fins, Cd ~ 0.65 (UW Madison)
    # ------------------------------------------------------------------------------

    # Converting to matlab file
    
        
    repository_root_path, _ = vehicle_parameters_functions.Get_Repository_Root_Path()    
    
    matlab_dict = {f"axial_array": axial_array, f"shear_array": shear_array, f"bending_array": bending_array, "length_along_rocket_linspace": length_along_rocket_linspace} # Dictionary to save as .mat file
    if location == "max_q":
        savemat(repository_root_path / "SFD" / "sfd_outputs_max_q.mat", matlab_dict) # Save as .mat file for MATLAB
    elif location == "off_the_rail":
        savemat(repository_root_path / "SFD" / "sfd_outputs_off_the_rail.mat", matlab_dict) # Save as .mat file for MATLAB
    # ------------------------------------------------------------------------------
    
    if print_inputs == True:
        print("Inputs:")
        if location == "max_q":
            print(f"\tWind gust at {location}: {wind_gust} m/s")
        elif location == "off_the_rail":
            print(f"\tEstimated rail whip at {location}: {wind_gust} m/s")
        print(f"\tAir density at {location}: {air_density} kg/m^3")
        print(f"\tAngle of attack at {location}: {AOA * (180 / np.pi):.2f} degrees")
        print(f"\tVelocity at {location}: {velocity} m/s")
        print(f"\tAcceleration at {location}: {acceleration} m/s^2")
        print(f"\tRocket mass at {location}: {total_mass:.2f} kg")
        print(f"\tCenter of c.GRAVITY at {location}: {CoM:.2f} m from bottom")
        print(f"\tCross sectional area: {S:.4f} m^2")
        print("-----------------------------------")

        print(f"Results at {location}:")
        print(f"\tMax shear force at {location}: {max(shear_array) * c.N2LBF:.2f} lbf")
        print(f"\tMax bending moment at {location}: {max(bending_array) * c.N2LBF * c.M2FT:.2f} lbf-ft")
        print(f"\tMax axial force at {location}: {max(axial_array) * c.N2LBF:.2f} lbf")
        print("-----------------------------------")

        print("Fin parameters:")
        print(f"\tFin center of pressure: {finCP:.2f} m from bottom")
        print(f"\tFin stability derivative: {finSD:.4f} ")
        print(f"\tFin lift: {finLift:.2f} N")
        print(f"\tRoot chord: {root_chord:.2f} m")
        print(f"\tTip chord: {tip_chord:.2f} m")
        print(f"\tSweep length: {sweep_length:.2f} m")
        print(f"\tFin height: {fin_height:.2f} m")
        print(f"\tNumber of fins: {numFins}")
        print("-----------------------------------")

        print("Nosecone parameters:")
        print(f"\tNose center of pressure: {noseCP:.2f} m from bottom")
        print(f"\tNose stability derivative: {noseSD:.4f} ")
        print(f"\tNose lift: {noseLift:.2f} N")
        print("-----------------------------------")

    if plot_on:
        # Plotting
        plt.figure()
        plot_num = 1
        for variable in ["shear_array", "bending_array", "axial_array"]:
            if variable == "shear_array":
                plot = shear_array * c.N2LBF
                ylabel = "Shear Force [lbf]"
                title = f"Shear Forces at {location}"
            if variable == "bending_array":
                plot = bending_array * c.N2LBF * c.M2FT
                ylabel = "Bending Moment [lbf-ft]"
                title = f"Bending Moments at {location}"
            if variable == "axial_array":
                plot = axial_array * c.N2LBF
                ylabel = "Axial Force [lbf]"
                title = f"Axial Forces at {location}"
            plt.subplot(1,3, plot_num)
            plt.plot(length_along_rocket_linspace * c.M2FT, plot)
            plt.title(title)
            plt.xlabel("Length from aft [ft]")
            plt.ylabel(ylabel)
            plot_num += 1
            plt.grid()

        plt.show()
            

def main():
    max_q_or_off_the_rail(plot_on = False, print_inputs = False)

if __name__ == "__main__":
    try:
        from SFD import parseWind
        from SFD import sfd
    except ModuleNotFoundError:
        import parseWind
        import sfd

    os.chdir(os.path.dirname(__file__))
    
    rerun_everything = False
    
    if rerun_everything:
        vehicle_main.vehicle_analysis()
    else:
        pass
        # parameters = csv_to_dataclass(parameters_csv_filepath)    
        
    max_q_or_off_the_rail(plot_on = True, print_inputs = True)
    