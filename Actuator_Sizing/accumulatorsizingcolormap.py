import main
import constants as c # type: ignore
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

volume_range = (10,100)
volume_step = 2
normalized_volume_range = ((int(volume_range[0]/volume_step), int(volume_range[1]/volume_step)))
cv_range = (0.5, 10)
cv_step = 0.1
fail_threshold = 0.25  # [s] actuation times above this (or stalls) are a fail
normalized_cv_range = ((int(cv_range[0]/cv_step), int(cv_range[1]/cv_step)))

actuation_times_grid = []

print(f"Total calcs to be performed: {(normalized_cv_range[1] - normalized_cv_range[0] + 1) * (normalized_volume_range[1] - normalized_volume_range[0] + 1)}")

for i in range(normalized_volume_range[0], normalized_volume_range[1] + 1, 1):
    volume = i * volume_step * c.IN2M**3
    actuation_times_line = []
    for j in range(normalized_cv_range[0], normalized_cv_range[1] + 1, 1):
        cv = j * cv_step
        Cv_total = main.sum_Cv(cv, main.Cv_solenoid)

        acc_volume_swept_history, acc_time_history, acc_angle_history, acc_time, acc_chamber_pressure_history, acc_flow_scfm_history, accumulator_pressure_history = main.actuation_time_kinematics_accumulator(
                main.rod_mass, main.piston_diameter, main.arm_length, main.friction_total, main.force_valve, Cv_total, main.pressure * c.PA2PSI, volume, main.dead_volume_m3, main.show_outputs, main.T_ambient_R, source_Cv=main.Cv_regulator)
        print(f"Actuation Time: {acc_time}s")
        actuation_times_line.append(acc_time)
    actuation_times_grid.append(actuation_times_line)

# Plotting the actuation times as a pass/fail heatmap
times = np.array(actuation_times_grid, dtype=float)  # rows = volume, columns = Cv
passing = np.isfinite(times) & (times <= fail_threshold)  # stalls come back as inf, so they count as fails

# Only passing cells keep a value. Masked (failing) cells are drawn in the colormap's "bad" colour, red.
pass_times = np.ma.masked_where(~passing, times)

# Gradient over the passing range only: bright green = fastest, dark green = slowest pass (right at the threshold)
pass_cmap = LinearSegmentedColormap.from_list("bright_to_dark_green", ["#00ff00", "#006400"])
pass_cmap.set_bad("red")
vmin = times[passing].min() if passing.any() else 0

fig, ax = plt.subplots()
im = ax.imshow(
    pass_times,
    # pad the extent by half a step so each cell is centred on its (Cv, volume) value
    extent=[cv_range[0] - cv_step / 2, cv_range[1] + cv_step / 2, volume_range[0] - volume_step / 2, volume_range[1] + volume_step / 2],
    origin='lower', aspect='auto', cmap=pass_cmap, vmin=vmin, vmax=fail_threshold, interpolation='nearest')
fig.colorbar(im, ax=ax, label=f'Actuation Time [s] (pass if <= {fail_threshold} s)')
ax.set_xticks(np.arange(cv_range[0], cv_range[1] + cv_step / 2, cv_step))
ax.set_yticks(np.arange(volume_range[0], volume_range[1] + volume_step / 2, volume_step))
ax.set_xlabel('Cv')
ax.set_ylabel('Accumulator Volume [in^3]')
ax.set_title('Actuation Time Heatmap')
ax.legend(handles=[Patch(facecolor='red', label=f'Fail (> {fail_threshold} s or stalled)')], loc='upper center', bbox_to_anchor=(0.5, -0.15))
fig.tight_layout()
plt.show()