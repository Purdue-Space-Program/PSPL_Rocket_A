import matplotlib.pyplot as plt
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
from matplotlib import animation
from dataclasses import dataclass

@dataclass
class MaterialClass:
    alpha: float

@dataclass
class NodeClass:
    i: int
    j: int
    temperature: float
    heat_transfer: float
    material: MaterialClass



def function_of_interest(multi_index):
    x = multi_index[1]
    y = multi_index[0]

    tank_wall_boundary_x_index = int((0.5 * x_width + 0.5)) # + 0.5 to go up to next integer if just below it, otherwise stays same integer

    if x < tank_wall_boundary_x_index:
        output = 10
    else:
        output = 100

    return(multi_index, output)



def Threaded_Run(function_of_interest, ndarray, USE_AI_SLOP):


    it = np.nditer(ndarray, flags=["multi_index"], op_flags=["readonly"])

    if USE_AI_SLOP:
        jobs = []
        for coordinate in it:
            jobs.append((it.multi_index, coordinate.copy()))

        with ThreadPoolExecutor() as executor:
            futures = [executor.submit(function_of_interest, idx) for idx, variable_input_combination in jobs]
            for f in tqdm(as_completed(futures), total=len(futures), desc="Threaded Run"):
                idx, output = f.result()

                ndarray[idx] = output

    else:
        for coordinate in tqdm(it, total=ndarray.size, desc="Not Threaded Run"):
            coordinate_index, output = function_of_interest(it.multi_index)
            ndarray[coordinate_index] = output

            # X, Y, Z = p.SetupArrays(variable_inputs_array, isp_map)
            # p.UpdateContinuousColorMap(X, Y, Z, constant_inputs_array)


def Setup_Arrays_for_Plotting(x_axis_name, y_axis_name, output_name, output_array):
    number_of_rows, number_of_columns = output_array.shape

    x = np.arange(number_of_columns)
    y = np.arange(number_of_rows)

    X, Y = np.meshgrid(x, y)

    return (X, Y, output_array)

def PlotColorMap(output_array, x_axis_name, y_axis_name, output_name, ax=None):
    if ax is None:
        ax = plt.gca()  # default to current axes

        X, Y, output_array = Setup_Arrays_for_Plotting(
                                                       x_axis_name,
                                                       y_axis_name,
                                                       output_name,
                                                       output_array,
                                                      )

    color_scheme = "RdBu_r"
    output_values_mesh = ax.pcolormesh(X, Y, output_array, cmap=color_scheme)
    # mesh = ax.contourf(X, Y, output_array, 100, cmap=color_scheme)


    # ax.set_title(f"{output_name.title()} of {inputs.constant_inputs['FUEL_NAME'][0].title()}", fontsize=8)
    ax.set_facecolor("lightgray")

    color_bar = plt.colorbar(output_values_mesh)
    # color_bar.set_label(output_label, fontsize=8)
    color_bar.ax.tick_params(labelsize=8)

    # # Make the colorbar use ~5 nice, round ticks over its value range
    color_bar.ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=3))
    color_bar.update_ticks()

    ax.set_xlabel("x", fontsize=8)
    ax.set_ylabel("y", fontsize=6)

    ax.tick_params(axis=x_axis_name, labelsize=8)  # X-axis tick numbers
    ax.tick_params(axis=y_axis_name, labelsize=8)  # Y-axis tick numbers

    # ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=6))
    # ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=6))



x_width, y_width = [3,10]
shape = y_width, x_width # fuck ah numpy indexing
ndarray = np.zeros(shape=shape)

Threaded_Run(function_of_interest, ndarray, False)

fig, ax = plt.subplots()
PlotColorMap(
             output_array=ndarray,
             x_axis_name="x",
             y_axis_name="y",
             output_name="output",
             ax=None,
            )
plt.show()

# t = np.linspace(0, 3, 40)
# g = -9.81
# v0 = 12
# z = g * t**2 / 2 + v0 * t

# v02 = 5
# z2 = g * t**2 / 2 + v02 * t

# scat = ax.scatter(t[0], z[0], c="b", s=5, label=f'v0 = {v0} m/s')
# line2 = ax.plot(t[0], z2[0], label=f'v0 = {v02} m/s')[0]
# ax.set(xlim=[0, 3], ylim=[-4, 10], xlabel='Time [s]', ylabel='Z [m]')
# ax.legend()


# def update(frame):
#     # for each frame, update the data stored on each artist.
#     x = t[:frame]
#     y = z[:frame]
#     # update the scatter plot:
#     data = np.stack([x, y]).T
#     scat.set_offsets(data)
#     # update the line plot:
#     line2.set_xdata(t[:frame])
#     line2.set_ydata(z2[:frame])
#     return (scat, line2)


# ani = animation.FuncAnimation(fig=fig, func=update, frames=40, interval=30)
# plt.show()