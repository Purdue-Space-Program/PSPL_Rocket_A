import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir(os.path.dirname(__file__))

'''
This script reads wind speed data from a text file, processes it to extract wind speed information,
and computes various statistics such as maximum, average, standard deviation, and percentiles of wind speed
measurements.

The file USW00094881-data.txt is a downloaded dataset from NOAA's National Centers for Environmental Information (NCEI).
The data was collected at Royal Center in Indiana, USA. The data was collected between 1918 and 1932.
This is the link to the dataset: https://www.ncei.noaa.gov/data/integrated-global-radiosonde-archive/access/data-por/
Search for USW00094881 to find the specific dataset.
'''

plot_on = True # Set to True to plot results, False to not plot

wind_data_file_path = "USW00094881-data.txt"
df = pd.read_csv(wind_data_file_path, sep='\s+', comment='#')
array_2d = df.values
# print(array_2d) # TEST
# print(array_2d[0:10]) # TEST
# print(array_2d[: , 8]) # Array of wind speeds TEST
wind_gust_speed_array = array_2d[: , 8] / 10 # Convert to m/s
max_wind_gust_speed = max(wind_gust_speed_array) # Max wind speed in m/s
avg_wind_gust_speed = np.mean(wind_gust_speed_array) # Average wind speed in m/s
std_dev_wind_gust_speed = np.std(wind_gust_speed_array) # Std dev of wind speed in m/s
percentile_75_wind_gust_speed = np.percentile(wind_gust_speed_array, 75) # 75th percentile wind speed in m/s
percentile_90_wind_gust_speed = np.percentile(wind_gust_speed_array, 90) # 90th percentile wind speed in m/s



if __name__ == "__main__":

    print(f"Max wind gust ever recorded: {max_wind_gust_speed:.2f} m/s") # Max wind speed in m/s
    print(f"Average wind gust: {avg_wind_gust_speed:.2f} m/s") # Average wind speed in m/s
    print(f"Standard deviation of wind gust speed: {std_dev_wind_gust_speed:.2f} m/s") # Std dev of wind speed in m/s
    print(f"75th percentile wind gust speed: {percentile_75_wind_gust_speed:.2f} m/s") # 75th percentile wind speed in m/s
    print(f"90th percentile wind gust speed: {percentile_90_wind_gust_speed:.2f} m/s") # 90th percentile wind speed in m/s

    if plot_on == True:
        plt.figure(figsize=(10, 6))
        plt.hist(wind_gust_speed_array, bins=20, alpha=0.7, edgecolor='black')
        plt.axvline(avg_wind_gust_speed, color='red', linestyle='--', linewidth=2, label=f'Mean: {avg_wind_gust_speed:.2f}')
        plt.axvline(percentile_75_wind_gust_speed, color='green', linestyle='--', linewidth=2, label=f'75th Percentile: {percentile_75_wind_gust_speed:.2f}')
        plt.axvline(percentile_90_wind_gust_speed, color='orange', linestyle='--', linewidth=2, label=f'90th Percentile: {percentile_90_wind_gust_speed:.2f}')
        plt.xlabel('Value')
        plt.ylabel('Frequency')
        plt.title('Distribution of Data')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()
