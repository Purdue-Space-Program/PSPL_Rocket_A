import matplotlib.pyplot as plt
import re
import numpy as np

starting_line_number = 16
ending_line_number = 56

actual_points_x = np.array([])
actual_points_y = np.array([])

goal_points_x = np.array([])
goal_points_y = np.array([])

actual_points_radius = np.array([])

with open("bulkhead_circularity_results.OUT", "r") as input_file:
    for current_line_number, current_line in enumerate(input_file, start=1):
        if (current_line_number >= starting_line_number) and (current_line_number <= ending_line_number):
             
            x_value_match = re.search(r"X([-+]?(?:\d+(?:\.\d*)?|\.\d+))", current_line)
            y_value_match = re.search(r"Y([-+]?(?:\d+(?:\.\d*)?|\.\d+))", current_line)

            if x_value_match and y_value_match:
                x_value = float(x_value_match.group(1))
                y_value = float(y_value_match.group(1))

                if current_line_number % 2 == 1:
                    actual_points_x = np.append(actual_points_x, x_value)
                    actual_points_y = np.append(actual_points_y, y_value)
                    
                    print(f"current_line_number: {current_line_number}")
                    print(f"current_line: {current_line}")
                    
                else:
                    goal_points_x = np.append(goal_points_x, x_value)
                    goal_points_y = np.append(goal_points_y, y_value)

                    # print(f"x_value: {x_value}")

goal_point_center_x = np.average(goal_points_x)
goal_point_center_y = np.average(goal_points_y)

actual_point_center_x = np.average(actual_points_x)
actual_point_center_y = np.average(actual_points_y)


actual_points_radius = (((actual_points_x - actual_point_center_x)**2) + ((actual_points_y - actual_point_center_y)**2))**(0.5)
actual_points_radius_average = np.average(actual_points_radius)
print(f"actual_points_radius_average: {actual_points_radius_average:.3f}")
print(f"actual_points_diameter_average: {2*actual_points_radius_average:.3f}")

plt.plot(actual_points_radius_average * np.cos(np.linspace(0, 6.28, num = 1000)) + goal_point_center_x, actual_points_radius_average * np.sin(np.linspace(0, 6.28, num = 1000)) + goal_point_center_y)

plt.scatter(actual_points_x, actual_points_y, color="red", zorder=2)
# plt.scatter(goal_points_x, goal_points_y, color="blue", zorder=1)
# plt.scatter(actual_points_x - goal_points_x, actual_points_y - goal_points_y, color="red", zorder=2)
plot_axes = plt.gca()
plot_axes.set_aspect("equal", adjustable="box")
plot_axes.set_box_aspect(1)
plt.grid()
plt.show()
