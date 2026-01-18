import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

'''
This script reads wind speed data from a text file, processes it to extract wind speed information,
and computes various statistics such as maximum, average, standard deviation, and percentiles of wind speed
measurements.

The file USW00094881-data.txt is a downloaded dataset from NOAA's National Centers for Environmental Information (NCEI).
The data was collected at Royal Center in Indiana, USA. The data was collected between 1918 and 1932.
This is the link to the dataset: https://www.ncei.noaa.gov/data/integrated-global-radiosonde-archive/access/data-por/
Search for USW00094881 to find the specific dataset.
'''

plot_on = False # Set to True to plot results, False to not plot

file_path  = 'USW00094881-data.txt'
df = pd.read_csv(file_path, sep='\s+', comment='#')
array_2d = df.values
print(array_2d)
print(array_2d[0:10])
print(array_2d[: , 8]) # Array of wind speeds
wspd_array = array_2d[: , 8] / 10 # Convert to m/s
max_wind_speed = max(wspd_array) # Max wind speed in m/s
avg_wind_speed = np.mean(wspd_array) # Average wind speed in m/s
std_dev_wind_speed = np.std(wspd_array) # Std dev of wind speed in m/s
percentile_75_wind_speed = np.percentile(wspd_array, 75) # 75th percentile wind speed in m/s
percentile_90_wind_speed = np.percentile(wspd_array, 90) # 90th percentile wind speed in m/s

print(f"Max wind speed ever recorded: {max_wind_speed} m/s") # Max wind speed in m/s
print(f"Average wind speed: {avg_wind_speed} m/s") # Average wind speed in m/s
print(f"Standard deviation of wind speed: {std_dev_wind_speed} m/s") # Std dev of wind speed in m/s
print(f"75th percentile wind speed: {percentile_75_wind_speed} m/s") # 75th percentile wind speed in m/s
print(f"90th percentile wind speed: {percentile_90_wind_speed} m/s") # 90th percentile wind speed in m/s

plt.figure(figsize=(10, 6))
plt.hist(wspd_array, bins=20, alpha=0.7, edgecolor='black')
plt.axvline(avg_wind_speed, color='red', linestyle='--', linewidth=2, label=f'Mean: {avg_wind_speed:.2f}')
plt.axvline(percentile_75_wind_speed, color='green', linestyle='--', linewidth=2, label=f'75th Percentile: {percentile_75_wind_speed:.2f}')
plt.axvline(percentile_90_wind_speed, color='orange', linestyle='--', linewidth=2, label=f'90th Percentile: {percentile_90_wind_speed:.2f}')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Distribution of Data')
plt.legend()
plt.grid(True, alpha=0.3)

if plot_on:
    plt.show()