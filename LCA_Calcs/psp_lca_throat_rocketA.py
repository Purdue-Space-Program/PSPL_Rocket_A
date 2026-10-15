
#python version of psp_lca_pathfinder_rocketA.m
#attempting to copy Chris's code for lca

#Basically the Biot number is still > 0.1 if we cut into hoops across throat section,
#so I also cut the thickness of the wall into hoops and now we have Biot<0.1

import os
import sys
import warnings

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.optimize import brentq

script_folder = os.path.dirname(os.path.abspath(__file__))
repo_folder   = os.path.dirname(script_folder)
data_folder   = os.path.join(script_folder, "Inputs")            # physical_change.csv
output_folder = os.path.join(script_folder, "Outputs_throat")
os.makedirs(output_folder, exist_ok=True)

sys.path.append(repo_folder)
import Vehicle_Level.vehicle_parameters as vehicle_parameters
parameters, _, _ = vehicle_parameters.main()


# Engine Contour
contour = np.loadtxt(os.path.join(repo_folder, "Chamber_Contour", "chamber_contour_meters.csv"), delimiter=",")
# Heat Transfer Coefficient
h_list = np.loadtxt(os.path.join(repo_folder, "Chamber_Temp", "chamber_heat_transfer_coefficient.csv"), delimiter=",")
# Gas Temp Along Wall
T_flame = np.loadtxt(os.path.join(repo_folder, "Chamber_Temp", "chamber_ambient_temp.csv"), delimiter=",")
# Pressure
wall_pressure = np.loadtxt(os.path.join(repo_folder, "Chamber_Temp", "chamber_pressure.csv"), delimiter=",")
# Temperature Derating of Material Properties
prop_file = os.path.join(data_folder, "physical_change.csv")
if os.path.isfile(prop_file):
    change_factor = np.loadtxt(prop_file, delimiter=",")
else:
    change_factor = None
    warnings.warn(f"{prop_file} not found. Yield strength is NOT derated with temperature "
                  "(room-temp value used) - stress results are unconservative.")

# MATERIAL PROPERTIES
# 1018 Steel (room temp values at T_REF)
STEEL_RHO = 7870            # [kg/m^3]      Density
STEEL_CP  = 486             # [J/(kg*K)]    Specific Heat
STEEL_K   = 51.9            # [W/(m*K)]     Thermal Conductivity
STEEL_S   = 2.4e8           # [Pa]          Yield Strength
STEEL_E   = 1.9e11          # [Pa]          Young's Modulus
STEEL_V   = 0.29            # [unitless]    Poisson's Ratio
STEEL_A   = 1.2e-5          # [1/K]         Thermal Expansion
T_REF     = 293             # [K]           Temp the room temp values are at

# Temperature dependence of each property
STEEL_K_SLOPE  = -3.24e-2                               # [W/(m*K^2)]   Conductivity slope (EN 1993-1-2, scaled to STEEL_K)
STEEL_A_SLOPE  = 0.8e-8                                 # [1/K^2]       Thermal expansion slope (EN 1993-1-2)
E_TEMPS        = np.array([294, 478, 589, 700, 811, 922])               # [K]
E_FACTOR       = np.array([1.0, 0.949, 0.901, 0.823, 0.694, 0.524])     # [unitless]    Young's modulus / room temp (ASME BPVC II-D TM-1, carbon steel)

# BURN PARAMETERS
T_BURN = parameters.burn_time   # [sec]     Burn Time
T_INIT = 300                    # [K]       Initial Temperature

# SAFETY PARAMETERS
T_MAX         = 773                     # [K]           Maximum Safe Temperature
SAFETY_FACTOR = 1.25                    # [unitless]    Thickness FoS
STRESS_FOS    = parameters.yield_FoS    # [unitless]    Yield FoS

# SLICING PARAMETERS
BI_LIMIT       = 0.1        # [unitless]    Max Biot number for a lump to be valid
N_SLICES_MIN   = 10         # [unitless]    Minimum slices through the wall of a sliced section
TEMP_CRITERION = "max"      # "max"  -> hottest point of wall (inner surface) <= T_MAX
                            # "mean" -> mean wall temp <= T_MAX (same criterion as plain LCA)
T_SEARCH       = (1e-4, 0.05)  # [m]        Thickness search range for sliced sections


# --- FUNCTIONS ---
def yield_strength(temp):
    # Yield derating from physical_change.csv, held flat outside its temp range
    if change_factor is None:
        return STEEL_S
    return STEEL_S * np.interp(temp, change_factor[:, 0], change_factor[:, 1])


def youngs_modulus(temp):
    return STEEL_E * np.interp(temp, E_TEMPS, E_FACTOR)


def thermal_expansion(temp):
    return STEEL_A + STEEL_A_SLOPE * (temp - T_REF)


def conductivity(temp):
    return STEEL_K + STEEL_K_SLOPE * (temp - T_REF)


def slice_wall(r_in, t_wall, h, T_gas):
    # Splits the wall into n annular slices, each thin enough to be lumped,
    # and returns the wall temperatures at the end of the burn.
    # Inner wall: convection from the gas. Outer wall: adiabatic (conservative).
    dr_max = BI_LIMIT * STEEL_K / h                         # [m] Slice thickness where slice Biot = BI_LIMIT
    n      = max(N_SLICES_MIN, int(np.ceil(t_wall / dr_max)))  # Number of slices
    dr     = t_wall / n                                     # [m] Slice thickness

    r_f = r_in + np.arange(n + 1) * dr                      # [m] Slice face radii
    r_c = (r_f[:-1] + r_f[1:]) / 2                          # [m] Slice centre radii

    # Per unit axial length
    C     = STEEL_RHO * STEEL_CP * np.pi * (r_f[1:]**2 - r_f[:-1]**2)               # [J/(m*K)] Slice heat capacity
    G     = 2*np.pi*STEEL_K / np.log(r_c[1:] / r_c[:-1])                            # [W/(m*K)] Slice-to-slice conductance
    G_gas = 1 / (1/(h*2*np.pi*r_in) + np.log(r_c[0]/r_in)/(2*np.pi*STEEL_K))        # [W/(m*K)] Gas-to-first-slice conductance

    K = np.diag(np.append(G, 0) + np.insert(G, 0, 0)) - np.diag(G, 1) - np.diag(G, -1)
    K[0, 0] += G_gas

    # C dT/dt = -K (T - T_gas)  ->  exact solution, no time stepping
    T = T_gas + expm(-(K / C[:, None]) * T_BURN) @ np.full(n, T_INIT - T_gas)

    return {
        "T_mean":   np.sum(C * T) / np.sum(C),                          # [K] Mean wall temp
        "T_hot":    T_gas - G_gas * (T_gas - T[0]) / (h*2*np.pi*r_in),  # [K] Inner surface temp
        "T_outer":  T[-1],                                              # [K] Outer wall temp
        "n":        n,
        "Bi_slice": h * dr / STEEL_K,
    }


def wall_temp(W):
    return W["T_mean"] if TEMP_CRITERION == "mean" else W["T_hot"]


def mark_hoop(ax, x_mid, hoop_index):
    if hoop_index is not None:
        ax.axvline(x_mid[hoop_index], color="black", linestyle=":")
        ax.text(x_mid[hoop_index], ax.get_ylim()[1], " Thick to Thin Wall Hoop Stress",
                rotation=90, va="top", fontsize=8)


def first_index(mask):
    idx = np.flatnonzero(mask)
    return idx[0] if idx.size else None


# --- CALCULATIONS ---

# Chamber Dimensions
x = contour[:, 0]      # [m] Depth
r = contour[:, 1]      # [m] Radius
N = len(x)             # Number of data points
M = N - 1              # Number of sections

# Surface Area
ds    = np.sqrt(np.diff(x)**2 + np.diff(r)**2)   # [m] Length of each segment
r_mid = (r[:-1] + r[1:]) / 2                     # [m] Mean radius of each segment
x_mid = (x[:-1] + x[1:]) / 2                     # [m] Midpoint of each segment
A_sec = 2*np.pi*r_mid * ds                       # [m^2] Surface area

# Sectioning Values
h_sec      = (h_list[:-1] + h_list[1:]) / 2      # [W/m^2-K] Sectioned h value
T_gas_wall = (T_flame[:-1] + T_flame[1:]) / 2    # [K] Section wall gas value

# - Single-Lump LCA (from chris)
log_term = np.log((T_gas_wall - T_INIT) / (T_gas_wall - T_MAX))

Lc_min   = h_sec * T_BURN / (STEEL_RHO * STEEL_CP * log_term)   # [m] Minimum Required Thickness
Lc_build = Lc_min * SAFETY_FACTOR                               # [m] Design Thickness

time_const = STEEL_RHO * STEEL_CP * Lc_build / h_sec
T_final    = T_gas_wall + (T_INIT - T_gas_wall) * np.exp(-T_BURN / time_const)

conduct = conductivity(T_final)

# Biot Number of the whole section
biot = h_sec * Lc_build / conduct

# - Sliced (Multi-Lump) Sections Where Single-Lump LCA Is Not Valid -
sliced     = biot > BI_LIMIT
feasible   = np.ones(M, dtype=bool)  # False -> no thickness keeps the wall under T_MAX
n_slices   = np.ones(M, dtype=int)   # [unitless] Slices through the wall
biot_slice = biot.copy()             # [unitless] Biot number of one slice
T_hot      = T_final.copy()          # [K] Inner (gas side) wall temp at end of burn
T_outer    = T_final.copy()          # [K] Outer wall temp at end of burn

for i in np.flatnonzero(sliced):
    wall = lambda t: slice_wall(r_mid[i], t, h_sec[i], T_gas_wall[i])
    over_limit = lambda t: wall_temp(wall(t)) - T_MAX

    if over_limit(T_SEARCH[1]) > 0:
        # Even a very thick wall gets too hot: heat-sink design fails here
        feasible[i] = False
        Lc_min[i]   = np.nan
        Lc_build[i] = np.nan
        W = wall(T_SEARCH[1])
    else:
        Lc_min[i]   = brentq(over_limit, *T_SEARCH, xtol=1e-7)
        Lc_build[i] = Lc_min[i] * SAFETY_FACTOR
        W = wall(Lc_build[i])

    T_final[i]    = W["T_mean"]
    T_hot[i]      = W["T_hot"]
    T_outer[i]    = W["T_outer"]
    n_slices[i]   = W["n"]
    biot_slice[i] = W["Bi_slice"]
    biot[i]       = h_sec[i] * Lc_build[i] / STEEL_K

# Recompute k at the hottest wall temp (same as T_final for single-lump sections)
conduct    = conductivity(T_hot)
biot_index = first_index(biot > BI_LIMIT)

# [kg] Engine Mass (annular volume, exact for thick walls at the throat)
sec_vol  = np.pi * ((r_mid + Lc_build)**2 - r_mid**2) * ds
sec_mass = STEEL_RHO * sec_vol

# -- HOOP STRESS --
# Inputs
r_i        = r_mid                                             # [m]   Inner radius
t_wall     = Lc_build                                          # [m]   Thickness
temp_dif   = T_gas_wall - T_final                              # [K]   Inner/Outer Wall Temp Diff
young_mod  = youngs_modulus(T_hot)                             # [Pa]  Young's Modulus
expansion  = thermal_expansion(T_hot)                          # [1/K] Thermal Expansion
max_stress = yield_strength(T_hot)                             # [Pa]  Maximum Allowable Stress

# Thermal stress: single-lump sections use the same estimate as the original
# script; sliced sections use the through-wall temp difference they solved for
q_therm      = h_sec * temp_dif
hoop_thermal = young_mod * expansion * q_therm * t_wall / (2 * (1 - STEEL_V) * conduct)
hoop_thermal[sliced] = (young_mod[sliced] * expansion[sliced]
                        * (T_hot[sliced] - T_outer[sliced]) / (2 * (1 - STEEL_V)))

# Pressure Calculations
pressure      = (wall_pressure[:-1] + wall_pressure[1:]) / 2
hoop_pressure = pressure * r_i / t_wall

# Validity
with np.errstate(invalid="ignore"):
    hoop_valid = r_i / t_wall
    hoop_index = first_index(hoop_valid < 10)

# Stress Calculations
rec_stress = max_stress / STRESS_FOS      # [Pa] Target stress
tot_stress = hoop_thermal + hoop_pressure

# --- PLOTTING ---
figs = {}

# Chamber Wall Thickness
fig, ax = plt.subplots(); figs["chamber_wall_thickness"] = fig
ax.plot(x_mid, Lc_build, linewidth=2, linestyle=":", color="red", label=f"Recommended FoS: {SAFETY_FACTOR}")
ax.plot(x_mid, Lc_min, linewidth=1.5, color="#d41e11", label="Absolute Minimum")
ax.axhline(np.nanmax(Lc_build), color="black", linestyle="--", label="Maximum Thickness")
if sliced.any():
    ax.axvspan(x_mid[sliced].min(), x_mid[sliced].max(), color="gray", alpha=0.15, label="Sliced (multi-lump)")
if (~feasible).any():
    ax.axvspan(x_mid[~feasible].min(), x_mid[~feasible].max(), color="red", alpha=0.15, label=f"No thickness meets {T_MAX} K")
ax.set_xlabel("Depth [m]"); ax.set_ylabel("Thickness [m]"); ax.set_title("Chamber Wall Thickness")
ax.grid(True); ax.legend(loc="upper left")

# Biot and Chamber Contour
fig, ax = plt.subplots(); figs["biot_number"] = fig
ax.plot(x_mid, biot, linewidth=2, label="Whole Wall (single lump)")
ax.plot(x_mid, biot_slice, linewidth=1.5, linestyle="--", label="One Slice (as used)")
if biot_index is not None:
    ax.axvline(x_mid[biot_index], linewidth=1.5, color="#4b0acc", label=f"Biot > {BI_LIMIT}")
ax.axhline(BI_LIMIT, color="red", linestyle="--", label="Max Valid Biot")
ax.set_xlabel("Chamber Depth [m]"); ax.set_ylabel("Biot Number [unitless]"); ax.set_title("Biot Number Along Contour")
ax.grid(True); ax.legend(loc="upper left")

# Built Chamber Section
fig, ax = plt.subplots(); figs["chamber_cross_section"] = fig
ax.plot(x, r, color="red", linestyle="--", label="Inner Contour")
ax.plot(x_mid, r_mid + Lc_build, color="black", linewidth=2, label="Outer Contour")
ax.plot(x, -r, color="red", linestyle="--")
ax.plot(x_mid, -r_mid - Lc_build, color="black", linewidth=2)
ax.set_aspect("equal")
ax.set_xlabel("Chamber Depth [m]"); ax.set_ylabel("Chamber Radius [m]"); ax.set_title("Chamber Cross Section")
ax.grid(True); ax.legend()

# Hoop, Thermal, Total Stress
for key, stress, label, title in [
    ("hoop_stress",    hoop_pressure, "Hoop Stress",    "Hoop Stress Over Chamber Contour"),
    ("thermal_stress", hoop_thermal,  "Thermal Stress", "Thermal Stress Over Chamber Contour"),
    ("total_stress",   tot_stress,    "Total Stress",   "Circumferential Stress Over Chamber Contour"),
]:
    fig, ax = plt.subplots(); figs[key] = fig
    ax.plot(x_mid, stress, color=(0.1, 0.4, 0.5) if key == "total_stress" else (0, 0.2, 0.9), linewidth=1, label=label)
    ax.plot(x_mid, max_stress, color=(0, 0.3, 0.7), linewidth=1, linestyle="--", label="Max Stress")
    ax.plot(x_mid, rec_stress, color=(0, 0.6, 0.9), linewidth=1, linestyle="--", label="Target Stress")
    mark_hoop(ax, x_mid, hoop_index)
    ax.set_xlabel("Chamber Depth [m]"); ax.set_ylabel("Stress [Pa]"); ax.set_title(title)
    ax.grid(True); ax.legend(loc="lower left")


# print
print(f"Sections with single-lump LCA: {np.sum(~sliced)}, sliced (multi-lump): {np.sum(sliced)}", end="")
if sliced.any():
    print(f" from x = {x_mid[sliced].min():.4f} to {x_mid[sliced].max():.4f} m")
    print(f"Slices per sliced section: {n_slices[sliced].min()} to {n_slices[sliced].max()}, "
          f"max slice Biot: {biot_slice[sliced].max():.3f} (limit {BI_LIMIT})")
    print(f"Thickness criterion in sliced sections: {TEMP_CRITERION} wall temp <= {T_MAX} K\n")
else:
    print("\n")

print(f"Max required thickness: {1000*np.nanmax(Lc_min):.2f} mm")
print(f"Max recommended thickness: {1000*np.nanmax(Lc_build):.2f} mm\n")

print(f"Max mean wall temp at end of burn with design thickness: {np.max(T_final[feasible]):.0f} K (limit {T_MAX} K)")
print(f"Max inner (gas side) wall temp at end of burn: {np.max(T_hot[feasible]):.0f} K")
print(f"Total wall mass at design thickness: {np.nansum(sec_mass):.2f} kg")

if (~feasible).any():
    print(f"\nWARNING: {np.sum(~feasible)} sections (x = {x_mid[~feasible].min():.4f} to "
          f"{x_mid[~feasible].max():.4f} m) cannot stay under {T_MAX} K at any thickness up to "
          f"{1000*T_SEARCH[1]:.0f} mm.\n         Inner wall reaches up to {T_hot[~feasible].max():.0f} K there. "
          "Thickness/stress/mass are left blank for these sections (mass total excludes them).")


#export
for name, fig in figs.items():
    fig.tight_layout()
    fig.savefig(os.path.join(output_folder, f"{name}.png"), dpi=200)

plt.show()
