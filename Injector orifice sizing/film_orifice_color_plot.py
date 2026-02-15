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

def calc_film_percent(area_total, dp):
    m_dot = area_total * cd * np.sqrt(2 * density * dp)
    film_percent = (m_dot / m_dot_ipa) * 100
    if film_percent > 30:
        return -1
    return film_percent

def display_color_plots():
    areas = np.linspace(0.001, 0.05, 2000) * IN22M2
    dp = np.linspace(1, 200, 2000) * PSI2PA

    X, Y = np.meshgrid(areas, dp)
    Z = np.zeros_like(X)
    print("This may take a while... please wait.")
    for i in range(len(areas)):
        for j in range(len(dp)):
            Z[i, j] = calc_film_percent(
                areas[j],
                dp[i]
            )
    Z[Z == -1] = np.nan

    percent_dp =  100 - (100 * chamber_pressure) / (chamber_pressure + Y)
    cmap = plt.cm.Reds.copy()
    cmap.set_bad('lightgray') # makes points with film % greater than 30 grey
    mesh = plt.pcolormesh(X * M22IN2, Y * PA2PSI, Z, cmap=cmap)
    cbar = plt.colorbar(mesh)
    cbar.set_label("Film %")
    plt.xlabel("Total Injection Area [in^2]")
    plt.ylabel("Pressure drop %")
    plt.title("Film %")
    plt.tight_layout()
    plt.show()

    cmap = plt.cm.Reds.copy()
    cmap.set_bad('lightgray') # makes points with film % greater than 30 grey
    mesh = plt.pcolormesh(X * M22IN2, percent_dp, Z, cmap=cmap)
    cbar = plt.colorbar(mesh)
    cbar.set_label("Film %")
    plt.xlabel("Total Injection Area [in^2]")
    plt.ylabel("Pressure drop %")
    plt.title(f"Film %\nChamber Pressure: {chamber_pressure * PA2PSI} psi, Cd = {cd}, Max Film = {max_film_percent}%")
    plt.tight_layout()
    plt.show()

display_color_plots()