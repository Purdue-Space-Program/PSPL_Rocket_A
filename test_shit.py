from CoolProp.CoolProp import PropsSI
import numpy as np
from pint import UnitRegistry
import constants as c
from dataclasses import dataclass
from itertools import product
import matplotlib.pyplot as plt
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

# u = UnitRegistry()

# nitrogen_standard_conditions_density = PropsSI("D", "P", 1 * c.ATM2PA, "T", c.T_AMBIENT, "nitrogen") * (u.kilogram / (u.meter**3))
# print(f"\nnitrogen_standard_conditions_density: {nitrogen_standard_conditions_density:.4g}\n")

# pressure_vessel_fluid_pressure = 500 * c.PSI2PA * u.pascal
# pressure_vessel_fluid_temperature = c.T_AMBIENT * u.kelvin
# pressure_vessel_fluid_density = PropsSI("D", "P", pressure_vessel_fluid_pressure.magnitude, "T", pressure_vessel_fluid_temperature.magnitude, "nitrogen") * (u.kilogram / (u.meter**3))
# pressure_vessel_orifice_area = np.pi * (((0.050 * c.IN2M) / 2)**2) * (u.meter**2)
# pressure_vessel_fluid_pressure_drop = pressure_vessel_fluid_pressure - 1 * c.ATM2PA * u.pascal

# def CalculateMassFlowRate(area, Cd, rho, pressure_drop):
#     m_dot = area * (Cd * np.sqrt(2 * rho * pressure_drop))
#     m_dot = m_dot.to(u.kilogram / u.second)
#     return m_dot


# pressure_vessel_fluid_mass_flow_rate = CalculateMassFlowRate(pressure_vessel_orifice_area,
#                                                              0.6,
#                                                              pressure_vessel_fluid_density,
#                                                              pressure_vessel_fluid_pressure_drop,
#                                                              )


# print(f"pressure_vessel_fluid_mass_flow_rate: {pressure_vessel_fluid_mass_flow_rate:.4g} kg/s")
# print(f"pressure_vessel_fluid_mass_flow_rate: {pressure_vessel_fluid_mass_flow_rate * c.KG2G:.4g} g/s\n")

# SCFM_2_KG_S = ((1**3 * c.FT32M3) * nitrogen_standard_conditions_density) / 60

# print(f"\nSCFM2KG_S: {SCFM_2_KG_S:.4g} kg/s")
# print(f"SCFM2KG_S: {SCFM_2_KG_S * c.KG2G:.4g} g/s\n")



#### CEA SHIT ######

def RunCEAWrap(
    chamber_pressure,
    OF_Ratio,
    fuel_name,
    oxidizer_name,
):

    # convert regular string for propellants to what CEA_wrap uses
    if fuel_name == "ethanol":
        CEA_fuel_name = cea_wrap.Fuel("C2H5OH(L)", temp=c.T_AMBIENT)
    elif fuel_name == "kerosene":
        CEA_fuel_name = cea_wrap.Fuel("Jet-A(L)", temp=c.T_AMBIENT)
    elif fuel_name == "ipa":
        CEA_fuel_name = cea_wrap.Fuel("C3H8O,2propanol", temp=c.T_AMBIENT)
    else:
        raise ValueError(f"{fuel_name} not supported")

    if oxidizer_name == "liquid oxygen":
        CEA_oxidizer_name = cea_wrap.Oxidizer("O2(L)", temp=90) # 90 K is temperature of oxidizer upon injection into combustion (same as copperhead's sizing)
    else:
        raise ValueError(f"{oxidizer_name} not supported")

    atmospheric_pressure = 1 * c.ATM2PA
    pressure_ratio = chamber_pressure / atmospheric_pressure # assume exit pressure is a constantly at the pressure of air a bit above sea level

    rocket = cea_wrap.RocketProblem(
        pressure =       chamber_pressure * c.PA2PSI,
        pip =            pressure_ratio, # pip is "Pressure ratio of chamber pressure to exit pressure." github.com/civilwargeeky/CEA_Wrap/blob/main/README.md#rocket-problem-constructor-additional-parameters
        materials =      [CEA_fuel_name, CEA_oxidizer_name],
        o_f =            OF_Ratio,
        pressure_units = "psi",
    )

    cea_results = rocket.run()
    
    return(cea_results)


def RunNASACEA(chamber_pressure,
                    OF_Ratio,
                    fuel_name,
                    oxidizer_name,
                ):
    
    # nasa's cea uses bar
    atmospheric_pressure = 1 * c.ATM2PA * c.PA2BAR
    chamber_pressure = chamber_pressure * c.PA2BAR


    # convert regular string for propellants to what CEA_wrap uses
    if fuel_name == "ethanol":
        CEA_fuel_name = "C2H5OH(L)"
        fuel_temperature = c.T_AMBIENT
    elif fuel_name == "kerosene":
        CEA_fuel_name = "Jet-A(L)"
        fuel_temperature = c.T_AMBIENT
    elif fuel_name == "ipa":
        CEA_fuel_name = "C3H8O,2propanol"
        fuel_temperature = c.T_AMBIENT
    else:
        raise ValueError(f"{fuel_name} not supported")

    if oxidizer_name == "liquid oxygen":
        CEA_oxidizer_name = "O2(L)"
        oxidizer_temperature = 90 # temperature of oxidizer upon injection into combustion (same as copperhead's sizing)
    else:
        raise ValueError(f"{oxidizer_name} not supported")


    reactant_names = [CEA_fuel_name, CEA_oxidizer_name]
    reactant_temperatures = np.array([fuel_temperature, oxidizer_temperature])

    fuel_weights = np.array([1, 0])
    oxidizer_weights = np.array([0, 1])

    reactants = nasa_cea.Mixture(reactant_names)
    products = nasa_cea.Mixture(reactant_names, products_from_reactants = True)

    reactant_weights = reactants.of_ratio_to_weights(
        oxidizer_weights,
        fuel_weights,
        OF_Ratio,
    )

    # needed to solve rocket problem
    chamber_enthalpy = (
        reactants.calc_property(
            nasa_cea.ENTHALPY,
            reactant_weights,
            reactant_temperatures,
        ) / nasa_cea.R) # divided by universal gas constant

    chamber_to_exit_pressure_ratio = chamber_pressure / atmospheric_pressure

    solver = nasa_cea.RocketSolver(products = products, reactants = reactants)
    solution = nasa_cea.RocketSolution(solver = solver)

    solver.solve(
        soln = solution,
        weights = reactant_weights,
        pc = chamber_pressure,
        pi_p = float(chamber_to_exit_pressure_ratio),
        hc = chamber_enthalpy,
        iac = True, # infinite area combustor
    )
    
    return(solution)

def SetupArrays(variable_inputs_array, x_axis_name, y_axis_name, output_name, output_array):
    x = np.array(variable_inputs_array[0, :][y_axis_name])
    y = np.array(variable_inputs_array[:, 0][x_axis_name])
    # z = np.array(output_array[output_name])

    Y, X = np.meshgrid(x, y) # I don't know why you have to swap X and Y but you do!
    # Z = z.reshape(len(x), len(y))
    # maybe change to this:
    # X, Y = np.meshgrid(x, y)
    # Z = z.reshape(len(y), len(x))

    Z = output_array
    return (X, Y, Z)


def PlotColorMap(variable_inputs_array, input_arrays, output_array, x_axis_name, y_axis_name, output_name, show_copv_limiting_factor, method_name, elapsed_time, ax=None):

    X, Y, output_values = SetupArrays(variable_inputs_array, x_axis_name, y_axis_name, output_name, output_array)

    if ax is None:
        ax = plt.gca()  # default to current axes

    axis_label_list, axis_values_list, output_label, output_values, color_scheme = FormatPlot(input_arrays, [x_axis_name, y_axis_name], [X, Y], output_name, output_values, show_copv_limiting_factor)

    X, Y = *axis_values_list,

    if show_copv_limiting_factor:
        output_values_mesh = ax.pcolormesh(X, Y, output_values, cmap=color_scheme, vmin=0, vmax=1*c.N2LBF)
    else:
        output_values_mesh = ax.pcolormesh(X, Y, output_values, cmap=color_scheme)
        # mesh = ax.contourf(X, Y, output_values, 100, cmap=color_scheme)


    # ax.set_title(f"{output_name.title()} of {inputs.constant_inputs['FUEL_NAME'][0].title()}", fontsize=8)
    ax.set_title(f"{method_name} elapsed time: {elapsed_time:.5g}", fontsize=12)
    ax.set_facecolor("lightgray")

    color_bar = plt.colorbar(output_values_mesh)
    color_bar.set_label(output_label, fontsize=12)
    color_bar.ax.tick_params(labelsize=12)

    # Make the colorbar use ~5 nice, round ticks over its value range
    color_bar.ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=3))
    color_bar.update_ticks()

    ax.set_xlabel(axis_label_list[0], fontsize=10)
    ax.set_ylabel(axis_label_list[1], fontsize=10)

    ax.tick_params(axis='x', labelsize=8)  # X-axis tick numbers
    ax.tick_params(axis='y', labelsize=8)  # Y-axis tick numbers

    ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=6))
    ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=6))


def FormatPlot(input_arrays, axis_name_list, axis_values_list, output_name, output_values, show_copv_limiting_factor):

    axis_label_list = [""] * len(axis_name_list)

    for input_array in input_arrays:
        for count, axis_name in enumerate(axis_name_list):
            if axis_name == input_array.name:
                axis_label_list[count] = input_array.plotting_axis_label
                axis_values_list[count] = axis_values_list[count] * input_array.plotting_axis_values_factor
            elif (input_array == input_arrays[-1]) and (axis_name == axis_name_list[-1]) and (axis_label_list[count] == ""):
                raise ValueError("axis name not recognized for plotting")
                
    
    # for input_array in enumerate(input_arrays):
    #     for count, axis_name in enumerate(axis_name_list):
    #         if axis_name == input_array.name:
    #             axis_label_list[count] = input_array.axis_label
    #             axis_values_list[count] = axis_values_list[count] * input_array.plotting_axis_values_factor
    #         elif (input_array == input_arrays[-1]) and (axis_name == axis_name_list[-1]):
    #             raise ValueError("axis name not recognized for plotting")
                
    
    # for count, axis_name in enumerate(axis_name_list):
    #     axis_values_factor = 1 # in case no factor is needed
    #     if axis_name == "OF_RATIO":
    #         axis_label = "OF Ratio"
    #     elif axis_name == "FUEL_TANK_LENGTH":
    #         axis_values_factor = c.M2IN
    #         axis_label = "Fuel Tank Length [in]"
    #     elif axis_name == "CONTRACTION_RATIO":
    #         axis_label = "Chamber to Throat Contraction Ratio"
    #     elif axis_name == "CHAMBER_PRESSURE":
    #         axis_values_factor = c.PA2PSI
    #         axis_label = "Chamber Pressure [psi]"
    #     else:
    #         raise ValueError("axis name not recognized for plotting")

    #     axis_label_list[count] = axis_label
    #     axis_values_list[count] = axis_values_list[count] * axis_values_factor

    # output_values_factor = 1 # in case no factor is needed
    # contour_lines = -1
    # if output_name == "JET_THRUST":
    #     output_values_factor = c.N2LBF
    #     output_label="Jet Thrust [lbf]"
    # elif output_name == "MASS_FLOW_RATE":
    #     output_label="Mass Flow Rate [kg/s]"
    # elif output_name == "ISP":
    #     # num_lines = 8
    #     # power = 1/4
    #     # contour_lines = np.max(output_array[output_name]) * 0.995 / (num_lines**power) * np.linspace(1, num_lines, num_lines)**power
    #     output_label="isp [s]"
    # elif output_name == "CHAMBER_TEMPERATURE":
    #     output_label="Chamber Temperature [k]"
    # elif output_name == "TANK_PRESSURE":
    #     output_values_factor = c.PA2PSI
    #     output_label="Tank Pressure [psi]"



    # elif output_name == "CHAMBER_INNER_DIAMETER":
    #     output_values_factor = c.M2IN
    #     output_label="Chamber Inner Diameter [in]"
    # elif output_name == "THROAT_DIAMETER":
    #     output_values_factor = c.M2IN
    #     output_label="Throat Diameter [in]"

    # elif output_name == "TOTAL_IMPULSE":
    #     output_label="Total Impulse [newtons-seconds]"
    # elif output_name == "APOGEE":
    #     output_values_factor = c.M2FT
    #     # altitude_limit = ax.contour(X, Y, Z, levels=[10000], colors='red', linewidths=2)
    #     # ax.clabel(altitude_limit, fmt='%d')
    #     output_label="Estimated Apogee [ft]"
    # elif output_name == "TAKEOFF_TWR":
    #     output_label="Takeoff TWR"
    # elif output_name == "RAIL_EXIT_TWR":
    #     output_label="Rail Exit TWR"
    # elif output_name == "WET_MASS":
    #     output_values_factor = c.KG2LB
    #     output_label="Wet Rocket Mass [lbm]"
    # elif output_name == "DRY_MASS":
    #     output_values_factor = c.KG2LB
    #     output_label="Dry Rocket Mass [lbm]"
    # elif output_name == "BURN_TIME":
    #     output_label="Burn Time [s]"
    # elif output_name == "RAIL_EXIT_VELOCITY":
    #     output_label="Rail Exit Velocity [m/s]"
    # elif output_name == "RAIL_EXIT_ACCELERATION":
    #     output_values_factor = 1 / c.GRAVITY
    #     output_label="Rail Exit Acceleration [G's]"
    # elif output_name == "MAX_ACCELERATION":
    #     output_values_factor = 1 / c.GRAVITY
    #     output_label="Max Acceleration [G's]"
    # elif output_name == "MAX_VELOCITY":
    #     output_values_factor = 1 / 343 # [m/s] 343 is da speed of sound
    #     output_label="Max Velocity [Mach]"

    # elif output_name == "TOTAL_LENGTH":
    #     output_values_factor = c.M2FT
    #     output_label="Total Length [ft]"

    # elif output_name == "OXIDIZER_TANK_LENGTH":
    #     output_values_factor = c.M2FT
    #     output_label="Oxidizer Tank Length [ft]"

    # elif output_name == "OXIDIZER_TANK_VOLUME":
    #     output_values_factor = c.M32L
    #     output_label="Oxidizer Tank Volume [liters]"

    # elif output_name == "OXIDIZER_TOTAL_MASS":
    #     output_values_factor = c.KG2LB
    #     output_label="Oxidizer Tank Mass [lbm]"

    # elif output_name == "FUEL_TANK_VOLUME":
    #     output_values_factor = c.M32L
    #     output_label="Fuel Tank Volume [liters]"

    # elif output_name == "FUEL_TOTAL_MASS":
    #     output_values_factor = c.KG2LB
    #     output_label="Fuel Tank Mass [lbm]"

    # elif output_name == "CHAMBER_STRAIGHT_WALL_LENGTH":
    #     output_values_factor = c.M2IN
    #     output_label="Chamber Straight Wall Length [in]"

    # elif output_name == "INJECTOR_TO_THROAT_LENGTH":
    #     output_values_factor = c.M2IN
    #     output_label="Injector To Throat Length [in]"

    # else:
    #     raise ValueError(f"{output_name} not recognized for plotting")


    # output_values = output_values * output_values_factor

    # # if contour_lines == -1:
    # #     ax.contour(X, Y, Z)
    # # else:
    # #     ax.contour(X, Y, Z, contour_lines)


    if show_copv_limiting_factor:
        color_scheme = "RdYlGn"

    else:
        color_scheme = "RdBu_r"
        "Spectral_r"

    output_label = "ISP [seconds]"
    return (axis_label_list, axis_values_list, output_label, output_values, color_scheme)




@dataclass
class InputArray:
    name: str
    array: np.ndarray
    datatype: type
    plotting_axis_values_factor: float
    plotting_unit: float
    
    def __post_init__(self):
        self.plotting_axis_label = f"{self.name} [{self.plotting_unit}]"

CEA_Wrap_unthreaded_index = 0
CEA_Wrap_threaded_index = 1
NASA_CEA_unthreaded_index = 2
NASA_CEA_threaded_index = 3

step_size_range = range(10, (50 + 1), 5)

output_arrays = np.full(shape = (4, len(step_size_range)),
                        fill_value = 0,
                        dtype = float,
                       )

number_of_tests = 10
for test_number in range(number_of_tests):
    
    for step_size_index, step_size in enumerate(step_size_range):

        chamber_pressure_array = InputArray(
                                            name = "Chamber pressure",
                                            array = np.linspace(100, 500, step_size) * c.PSI2PA,
                                            datatype = np.float32,
                                            plotting_axis_values_factor = c.PA2PSI,
                                            plotting_unit = "PSI"
                                            )
        OF_ratio_array = InputArray(
                                            name = "OF ratio",
                                            array = np.linspace(0.5, 2.5, step_size),
                                            datatype = np.float32,
                                            plotting_axis_values_factor = 1,
                                            plotting_unit = "N/A"
                                            )
        # contraction_ratio_array = InputArray(
        #                                     name = "Contraction ratio",
        #                                     array = np.linspace(2, 15, step_size),
        #                                     datatype = np.float32,
        #                                     plotting_axis_values_factor = 1,
        #                                     plotting_unit = "N/A"
        #                                     )

        # input_arrays = (chamber_pressure_array, OF_ratio_array, contraction_ratio_array)
        input_arrays = (OF_ratio_array, chamber_pressure_array)

        shape = [step_size]*len(input_arrays)
        fields_dtype = []
        for input_array in input_arrays:
            fields_dtype.append((input_array.name, input_array.datatype))


        possible_combinations = list(product(*(input_array.array for input_array in input_arrays)))
        # holy shit i cooked
        possible_combinations = np.array(possible_combinations, dtype=np.dtype(fields_dtype))
        possible_combinations = possible_combinations.reshape(shape)

        # lowercase strings here because you can't in an ndarray



        # use_NASA_cea = True


        import CEA_Wrap as cea_wrap
        import cea as nasa_cea

        def NASACEAJob(current_index, combination):
            
            cea_results = RunNASACEA(
                            chamber_pressure = combination["Chamber pressure"],
                            OF_Ratio = combination["OF ratio"],
                            fuel_name = "ipa", 
                            oxidizer_name = "liquid oxygen"
                            )

            return(current_index, cea_results)

        def CEAWrapJob(current_index, combination):
            
            cea_results = RunCEAWrap(
                            chamber_pressure = combination["Chamber pressure"],
                            OF_Ratio = combination["OF ratio"],
                            fuel_name = "ipa", 
                            oxidizer_name = "liquid oxygen"
                            )

            return(current_index, cea_results)



        def RunCEA(use_NASA_CEA, use_threading):
            
            elapsed_time = np.nan # to establish its existence out of the tqdm context manager so it doesn't delete after i assign it within the tqdm context manager
            
            output_array = np.full(
                shape=np.shape(possible_combinations),
                fill_value=np.nan,
                dtype=float,
            )

            if use_NASA_CEA == True:
                description = "NASA CEA"
                
                
                if use_threading == False:
                    with tqdm(np.ndenumerate(possible_combinations), total=possible_combinations.size, desc=f"{description} unthreaded run") as progress_bar: # np.ndenumerate is a bit faster than using np.nditer, and also a bit cleaner
                        for current_index, combination in progress_bar: 
                            cea_results = RunNASACEA(
                                                chamber_pressure = combination["Chamber pressure"],
                                                OF_Ratio = combination["OF ratio"],
                                                fuel_name = "ipa",
                                                oxidizer_name = "liquid oxygen",
                                                )
                                    
                            # [-1] give the last station which is the exit nozzle
                            output_isp = cea_results.Isp[-1] / c.STANDARD_GRAVITY
                            # exit_pressure = cea_results.P[-1]

                            output_array[current_index] = output_isp
                            
                        elapsed_time = progress_bar.format_dict["elapsed"]

                elif use_threading == True:
                    jobs = []
                    
                    # for current_index, combination in np.ndenumerate(possible_combinations):
                    possible_combinations_iterator = np.nditer(possible_combinations, flags=["multi_index"], op_flags=["readonly"],) # need this for tqdm to work
                    
                    for combination in tqdm(possible_combinations_iterator, total=possible_combinations.size, desc="Create Jobs",):
                        current_index = possible_combinations_iterator.multi_index
                    
                        jobs.append((current_index, combination.copy()))

                    with ThreadPoolExecutor() as executor:
                        futures = [executor.submit(NASACEAJob, current_index, combination) for current_index, combination in jobs]
                        
                        
                        with tqdm(as_completed(futures), total=possible_combinations.size, desc=f"{description} threaded run") as progress_bar:
                            for future in progress_bar:
                                current_index, cea_results = future.result()
                                output_array[current_index] = cea_results.Isp[-1] / c.STANDARD_GRAVITY
                        
                        elapsed_time = progress_bar.format_dict["elapsed"]

            elif use_NASA_CEA == False:
                description = "CEA Wrap"

                if use_threading == False:
                    
                    with tqdm(np.ndenumerate(possible_combinations), total=possible_combinations.size, desc=f"{description} unthreaded run") as progress_bar: # np.ndenumerate is a bit faster than using np.nditer, and also a bit cleaner
                        for current_index, combination in progress_bar: 
                            
                            cea_results = RunCEAWrap(
                                                chamber_pressure = combination["Chamber pressure"],
                                                OF_Ratio = combination["OF ratio"],
                                                fuel_name = "ipa", 
                                                oxidizer_name = "liquid oxygen",
                                                )
                            output_array[current_index] = cea_results.isp
                        
                        elapsed_time = progress_bar.format_dict["elapsed"]
                    
                elif use_threading == True:
                    jobs = []
                    
                    # for current_index, combination in np.ndenumerate(possible_combinations):
                    possible_combinations_iterator = np.nditer(possible_combinations, flags=["multi_index"], op_flags=["readonly"],) # need this for tqdm to work
                    
                    for combination in tqdm(possible_combinations_iterator, total=possible_combinations.size, desc="Create Jobs",):
                        current_index = possible_combinations_iterator.multi_index
                    
                        jobs.append((current_index, combination.copy()))

                    with ThreadPoolExecutor() as executor:
                        futures = [executor.submit(CEAWrapJob, current_index, combination) for current_index, combination in jobs]
                        
                        
                        with tqdm(as_completed(futures), total=possible_combinations.size, desc=f"{description} threaded run") as progress_bar:
                            for future in progress_bar:
                                current_index, cea_results = future.result()
                                output_array[current_index] = cea_results.isp
                            
                        elapsed_time = progress_bar.format_dict["elapsed"]
                    
                    # print(f"elapsed_time: {elapsed_time}")

            return(output_array, elapsed_time)

        # CEA_Wrap_unthreaded_output_array, CEA_Wrap_unthreaded_elapsed_time = RunCEA(use_NASA_CEA = False, use_threading = False)
        # CEA_Wrap_threaded_output_array, CEA_Wrap_threaded_elapsed_time = RunCEA(use_NASA_CEA = False, use_threading = True)
        NASA_CEA_unthreaded_output_array, NASA_CEA_unthreaded_elapsed_time = RunCEA(use_NASA_CEA = True, use_threading = False)
        NASA_CEA_threaded_output_array, NASA_CEA_threaded_elapsed_time = RunCEA(use_NASA_CEA = True, use_threading = True)


        # output_arrays[CEA_Wrap_unthreaded_index].append(CEA_Wrap_unthreaded_elapsed_time)
        # output_arrays[CEA_Wrap_threaded_index].append(CEA_Wrap_threaded_elapsed_time)
        output_arrays[(NASA_CEA_unthreaded_index, step_size_index)] += NASA_CEA_unthreaded_elapsed_time
        output_arrays[(NASA_CEA_threaded_index, step_size_index)] += NASA_CEA_threaded_elapsed_time

        # fig, axs = plt.subplots(2, 2, figsize=(18,10), constrained_layout=False)
        # fig.suptitle(f"CEA_Wrap and NASA CEA speed comparison", fontsize=14)
        # PlotColorMap(possible_combinations, input_arrays, CEA_Wrap_unthreaded_output_array, "OF ratio", "Chamber pressure", "ISP", False, method_name = "CEA_Wrap_unthreaded", elapsed_time = CEA_Wrap_unthreaded_elapsed_time, ax=axs[0, 0])
        # PlotColorMap(possible_combinations, input_arrays, CEA_Wrap_threaded_output_array, "OF ratio", "Chamber pressure", "ISP", False, method_name = "CEA_Wrap_threaded", elapsed_time = CEA_Wrap_threaded_elapsed_time, ax=axs[1, 0])
        # PlotColorMap(possible_combinations, input_arrays, NASA_CEA_unthreaded_output_array, "OF ratio", "Chamber pressure", "ISP", False, method_name = "NASA_CEA_unthreaded", elapsed_time = NASA_CEA_unthreaded_elapsed_time, ax=axs[0, 1])
        # PlotColorMap(possible_combinations, input_arrays, NASA_CEA_threaded_output_array, "OF ratio", "Chamber pressure", "ISP", False, method_name = "NASA_CEA_threaded", elapsed_time = NASA_CEA_threaded_elapsed_time, ax=axs[1, 1])

        # plt.subplots_adjust(hspace = 0.3, left = 0.08, right = 0.96)
        # plt.show()

        # run_timed_test = True

        # if run_timed_test:
        #     import timeit
        #     timed_tests = 3
        #     CEA_Wrap_elapsed_time = timeit.timeit(CEAWrapLoop, number=timed_tests)
        #     NASA_CEA_elapsed_time = timeit.timeit(NASACEALoop, number=timed_tests)
        #     print(f"Number of timed tests: {timed_tests}")
        #     print(f"Average CEA_Wrap_elapsed_time: {CEA_Wrap_elapsed_time/timed_tests}")
        #     print(f"Average NASA_CEA_elapsed_time: {NASA_CEA_elapsed_time/timed_tests}")

output_arrays = output_arrays / number_of_tests

elapsed_time_figure, elapsed_time_axes = plt.subplots(1, 1, figsize=(18,10), constrained_layout=False)
elapsed_time_axes.set_title(f"Number of tests: {number_of_tests}")
# elapsed_time_axes.plot((np.asarray(step_size_range))**2, output_arrays[CEA_Wrap_unthreaded_index], label = "CEA Wrap unthreaded")
# elapsed_time_axes.plot((np.asarray(step_size_range))**2, output_arrays[CEA_Wrap_threaded_index], label = "CEA Wrap threaded")
elapsed_time_axes.plot((np.asarray(step_size_range))**2, output_arrays[NASA_CEA_unthreaded_index], label = "NASA CEA unthreaded")
elapsed_time_axes.plot((np.asarray(step_size_range))**2, output_arrays[NASA_CEA_threaded_index], label = "NASA CEA threaded")
elapsed_time_axes.legend()
elapsed_time_axes.set_ylabel("Time [seconds]")
elapsed_time_axes.set_xlabel("# of iterations")
plt.show()
