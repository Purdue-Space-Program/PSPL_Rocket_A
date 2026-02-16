import numpy as np
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import vehicle_parameters as v
import matplotlib.pyplot as plt
from constants import *

m_dot_ipa = v.parameters.core_fuel_mass_flow_rate
density = DENSITY_IPA
cd = 0.6
max_film_percent = 30
max_m_dot = m_dot_ipa * max_film_percent / 100
chamber_pressure = v.parameters.chamber_pressure 

def display_color_plots():
    areas = np.linspace(0.005, 0.03, 2000) * IN22M2
    dp = np.linspace(chamber_pressure / 0.8 - chamber_pressure, 200 * PSI2PA, 2000)

    X, Y = np.meshgrid(areas, dp)
    Z = np.zeros_like(X)
    print("This may take a while... please wait.")
    m_dot = X * cd * np.sqrt(2 * density * Y)
    Z = (m_dot / m_dot_ipa) * 100
    Z[(Z > 30) | (Z < 15)] = np.nan

    valid_mask = (Z >= 15) & (Z <= 30)
    
    film_ranges = []
    for col in range(Z.shape[1]):  # Iterate over each area
        valid_values = Z[:, col][valid_mask[:, col]]
        if len(valid_values) > 0:
            film_range = np.max(valid_values) - np.min(valid_values)
            min_film = np.min(valid_values)
            max_film = np.max(valid_values)
        else:
            film_range = 0
            min_film = np.nan
            max_film = np.nan
        film_ranges.append((film_range, min_film, max_film))
    
    film_ranges = np.array(film_ranges)
    
    # Find the area with the largest range
    best_col = np.argmax(film_ranges[:, 0])
    best_area = areas[best_col] * M22IN2
    best_range = film_ranges[best_col, 0]
    best_min_film = film_ranges[best_col, 1]
    best_max_film = film_ranges[best_col, 2]
    col_data = Z[:, best_col]
    valid_col = valid_mask[:, best_col]
    
    idx_min = np.where((col_data == best_min_film) & valid_col)[0][0]
    idx_max = np.where((col_data == best_max_film) & valid_col)[0][0]
    
    dp_min = dp[idx_min] * PA2PSI
    dp_max = dp[idx_max] * PA2PSI


    percent_dp =  100 - (100 * chamber_pressure) / (chamber_pressure + Y)
    percent_dp_min = percent_dp[idx_min, best_col]
    percent_dp_max = percent_dp[idx_max, best_col]
    plt.axvline(best_area, color='yellow', linestyle='--', linewidth=2, alpha=0.7,
                label=f'Best Area: {best_area:.4f} in² (Range: {best_range:.1f}%)')
    plt.plot(best_area, dp_min, 'go', markersize=10, markeredgecolor='black', 
             label=f'Min: {best_min_film:.1f}% @ {dp_min:.1f} psi')
    plt.plot(best_area, dp_max, 'ro', markersize=10, markeredgecolor='black',
             label=f'Max: {best_max_film:.1f}% @ {dp_max:.1f} psi')
    
    cmap = plt.cm.coolwarm.copy()
    cmap.set_bad('lightgray') # makes points with film % greater than 30 grey
    mesh = plt.pcolormesh(X * M22IN2, Y * PA2PSI, Z, cmap=cmap)
    cbar = plt.colorbar(mesh)
    cbar.set_label("Film %")
    plt.xlabel("Total Injection Area [in^2]")
    plt.ylabel("Pressure drop [psi]")
    plt.title("Film %")
    plt.tight_layout()
    plt.legend()
    plt.show()

    cmap = plt.cm.coolwarm.copy()
    cmap.set_bad('lightgray') # makes points with film % greater than 30 grey
    mesh = plt.pcolormesh(X * M22IN2, percent_dp, Z, cmap=cmap)
    cbar = plt.colorbar(mesh)
    cbar.set_label("Film %")

    plt.axvline(best_area, color='yellow', linestyle='--', linewidth=2, alpha=0.7,
                label=f'Best Area: {best_area:.4f} in² (Range: {best_range:.1f}%)')
    plt.plot(best_area, percent_dp_min, 'go', markersize=10, markeredgecolor='black',
             label=f'Min: {best_min_film:.1f}% @ {percent_dp_min:.1f}%')
    plt.plot(best_area, percent_dp_max, 'ro', markersize=10, markeredgecolor='black',
             label=f'Max: {best_max_film:.1f}% @ {percent_dp_max:.1f}%')
    
    plt.xlabel("Total Injection Area [in^2]")
    plt.ylabel("Pressure drop %")
    plt.title(f"Film %\nChamber Pressure: {chamber_pressure * PA2PSI} psi, Cd = {cd}, Max Film = {max_film_percent}%")
    plt.legend()
    plt.tight_layout()
    plt.show()

display_color_plots()