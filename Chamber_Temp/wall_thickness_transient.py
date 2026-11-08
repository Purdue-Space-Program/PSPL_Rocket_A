#Trial wall thickness from the 400 K isotherm, then a 1D transient conduction check on that thickness
#(no LCA here, the temperature varies through the wall)
#
#Step 1: run the Bartz wall conduction solve on the current design wall, find how deep the 400 K
#        isotherm gets at the end of the burn -> that depth is the trial wall thickness (saved to csv)
#        the thicker wall is the trial wall with a manufacturable straight section, blended back on before the throat
#Step 2: re-solve the temperature through the thicker wall with the 1D transient heat equation,
#        infinite cylinder of radius r_o = L (the wall thickness):
#           hot face (r* = 1):  convection from the gas, h (T_aw - T)
#           outer face (r* = 0): dT/dr = 0 (adiabatic, heat stays in the wall, conservative)
#        exact series solution (Incropera 5.47):
#           theta* = (T - T_aw) / (T_i - T_aw) = sum C_n exp(-zeta_n^2 Fo) J0(zeta_n r*)
#           zeta_n J1(zeta_n) / J0(zeta_n) = Bi,  C_n = (2 / zeta_n) J1(zeta_n) / (J0^2(zeta_n) + J1^2(zeta_n))
#           Bi = h L / k,  Fo = alpha t / L^2,  r* = (L - depth) / L
#        
#Step 3: thermal + pressure stresses in the thicker wall
#        pressure hoop stress both ways: thin wall (p r / t) and thick wall (Lame),
#        thin wall is only valid where r / t >= THIN_WALL_RATIO

import os
import sys
import importlib.util

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.special import j0, j1, jn_zeros

script_folder = os.path.dirname(os.path.abspath(__file__))
repo_folder   = os.path.dirname(script_folder)
output_folder = os.path.join(script_folder, "Outputs_transient")
os.makedirs(output_folder, exist_ok=True)

sys.path.append(repo_folder)
import Vehicle_Level.vehicle_parameters as vehicle_parameters
import constants as c
parameters, _, _ = vehicle_parameters.main()

#reusing the wall conduction solver + isotherm finder from the Bartz script so the isotherm matches its plot
bartz_spec = importlib.util.spec_from_file_location("bartz", os.path.join(script_folder, "Combustion Chamber Bartz Equation.py"))
bartz = importlib.util.module_from_spec(bartz_spec)
bartz_spec.loader.exec_module(bartz)

# INPUTS (from the Bartz script csvs)
contour       = np.loadtxt(os.path.join(repo_folder, "Chamber_Contour", "chamber_contour_meters.csv"), delimiter=",")
h_gas         = np.loadtxt(os.path.join(script_folder, "chamber_heat_transfer_coefficient.csv"), delimiter=",")  # [W/m^2-K]
T_aw          = np.loadtxt(os.path.join(script_folder, "chamber_ambient_temp.csv"), delimiter=",")               # [K] Adiabatic wall temp
wall_pressure = np.loadtxt(os.path.join(script_folder, "chamber_pressure.csv"), delimiter=",")                   # [Pa]

prop_file = os.path.join(repo_folder, "LCA_Calcs", "Inputs", "physical_change.csv")
change_factor = np.loadtxt(prop_file, delimiter=",")   # [K, yield / room temp yield]

# MATERIAL PROPERTIES
# 1018 Steel (room temp values)
STEEL_RHO = 7870            # [kg/m^3]      Density
STEEL_CP  = 486             # [J/(kg*K)]    Specific Heat
STEEL_K   = 51.9            # [W/(m*K)]     Thermal Conductivity
STEEL_S   = 2.4e8           # [Pa]          Yield Strength
STEEL_E   = 1.9e11          # [Pa]          Young's Modulus
STEEL_V   = 0.29            # [unitless]    Poisson's Ratio
STEEL_A   = 1.2e-5          # [1/K]         Thermal Expansion
T_REF     = 293             # [K]           Temp the room temp values are at
STEEL_ALPHA = STEEL_K / (STEEL_RHO * STEEL_CP)  # [m^2/s] Thermal diffusivity

STEEL_A_SLOPE = 0.8e-8                                              # [1/K^2]   Thermal expansion slope (EN 1993-1-2)
E_TEMPS       = np.array([294, 478, 589, 700, 811, 922])            # [K]
E_FACTOR      = np.array([1.0, 0.949, 0.901, 0.823, 0.694, 0.524])  # [unitless] Young's modulus / room temp (ASME BPVC II-D TM-1)

# PARAMETERS
ISOTHERM_TEMP  = 400                    # [K]       Isotherm that sets the trial thickness
THICKNESS_FOS  = 1.0                    # [unitless] Multiplier on the isotherm depth (1.0 = thickness is exactly the isotherm depth)
MIN_THICKNESS  = 1.5e-3                 # [m]       Floor where the isotherm is shallow/absent (nozzle exit)
STRAIGHT_WALL_THICKNESS = 0.25 * c.IN2M # [m]       Manufacturable thickness for the straight wall section
BLEND_LENGTH   = 0.5 * c.IN2M           # [m]       Axial length of the blend from the straight wall thickness onto the trial thickness
T_BURN         = parameters.burn_time   # [s]       Burn time
T_INIT         = 300                    # [K]       Initial wall temp (same as the Bartz script)
STRESS_FOS     = parameters.yield_FoS   # [unitless] Yield FoS
N_TERMS        = 100                    # Series terms
N_DEPTH        = 80                     # Points through the wall
THIN_WALL_RATIO = 10                    # [unitless] r / t above which thin wall hoop stress is valid



def yield_strength(temp):
    return STEEL_S * np.interp(temp, change_factor[:, 0], change_factor[:, 1])


def youngs_modulus(temp):
    return STEEL_E * np.interp(temp, E_TEMPS, E_FACTOR)


def thermal_expansion(temp):
    return STEEL_A + STEEL_A_SLOPE * (temp - T_REF)


J1_ZEROS = np.concatenate(([0], jn_zeros(1, N_TERMS - 1)))   # roots of J1, starting at 0
J0_ZEROS = jn_zeros(0, N_TERMS)                               # roots of J0


def eigenvalues(Bi, n_terms=N_TERMS):
    # roots of zeta J1(zeta) / J0(zeta) = Bi, the n-th root sits between the (n-1)-th root of J1 and the n-th root of J0
    f = lambda z: z * j1(z) - Bi * j0(z)
    return np.array([brentq(f, J1_ZEROS[n], J0_ZEROS[n]) for n in range(n_terms)])


def wall_transient(h, T_gas, L, depths, t):
    # 1D transient infinite cylinder of radius L, convection on the hot face (r* = 1), adiabatic outer face (r* = 0)
    # depths: [m] distance into the wall from the hot face, t: [s] time (scalar or array)
    Bi   = h * L / STEEL_K
    zeta = eigenvalues(Bi)
    C_n  = 2 / zeta * j1(zeta) / (j0(zeta)**2 + j1(zeta)**2)
    Fo   = STEEL_ALPHA * np.atleast_1d(t) / L**2

    r_star = (L - depths) / L                                   # 0 at the adiabatic face, 1 at the hot face
    theta  = (C_n[None, None, :] * np.exp(-zeta[None, None, :]**2 * Fo[:, None, None])
              * j0(zeta[None, None, :] * r_star[None, :, None])).sum(axis=2)
    T = T_gas + (T_INIT - T_gas) * theta                        # [K] shape (n_times, n_depths)
    return T if np.ndim(t) else T[0]


def cylinder_stresses(r, T, p, E, a_exp):
    # thick walled cylinder, inner radius r[0], outer radius r[-1]
    # thermal: Timoshenko & Goodier, generalized plane strain (free ends)
    # pressure: Lame with closed ends
    a, b = r[0], r[-1]
    Tr        = T * r
    I_r       = np.concatenate(([0], np.cumsum((Tr[1:] + Tr[:-1]) / 2 * np.diff(r))))  # int_a^r T r dr
    I_b       = I_r[-1]
    K_th      = E * a_exp / (1 - STEEL_V)
    with np.errstate(invalid="ignore", divide="ignore"):
        sig_r_th = K_th / r**2 * ((r**2 - a**2) / (b**2 - a**2) * I_b - I_r)
        sig_t_th = K_th / r**2 * ((r**2 + a**2) / (b**2 - a**2) * I_b + I_r - T * r**2)
    sig_z_th  = K_th * (2 * I_b / (b**2 - a**2) - T)

    lame      = p * a**2 / (b**2 - a**2)
    sig_r_p   = lame * (1 - b**2 / r**2)
    sig_t_p   = lame * (1 + b**2 / r**2)
    sig_z_p   = np.full_like(r, lame)

    s_r, s_t, s_z = sig_r_th + sig_r_p, sig_t_th + sig_t_p, sig_z_th + sig_z_p
    von_mises = np.sqrt(((s_r - s_t)**2 + (s_t - s_z)**2 + (s_z - s_r)**2) / 2)
    return {"t_th": sig_t_th, "t_p": sig_t_p, "t_tot": s_t, "vm": von_mises}


# --- CALCULATIONS ---
x_wall   = contour[:, 0]                 # [m] Axial position
r_inner  = contour[:, 1]                 # [m] Inner radius
throat   = np.argmin(r_inner)
x_rel    = x_wall - x_wall[throat]       # [m] Axial position relative to throat
N        = len(x_wall)

# 400 K isotherm depth on the current design wall
design_wall = (parameters.chamber_outer_diameter - parameters.chamber_inner_diameter) / 2  # [m]
fv_depths, fv_temps = bartz.wall_temperature_field(h_gas, T_aw, r_inner, design_wall, STEEL_K, T_BURN)
iso_depth = bartz.isotherm_depth(fv_depths, fv_temps, ISOTHERM_TEMP)

capped = iso_depth >= design_wall
trial_thickness = np.maximum(iso_depth * THICKNESS_FOS, MIN_THICKNESS)   # [m]

trial_outer_diameter = 2 * (r_inner + trial_thickness)   # [m] Outer diameter of the trial wall

# Thicker wall: isotherm depth is too thin to manufacture in the straight wall, so hold STRAIGHT_WALL_THICKNESS
# from the injector until the isotherm catches up, then a cubic blend (matches value + slope at both ends) onto the trial thickness.
# The throat and everything downstream stay at the trial thickness
catch_up     = np.argmax((np.arange(N) < throat) & (trial_thickness >= STRAIGHT_WALL_THICKNESS))  # first station upstream of the throat at the target
blend_start  = x_wall[catch_up] - BLEND_LENGTH / 2                                                # [m]
blend_end_i  = np.searchsorted(x_wall, x_wall[catch_up] + BLEND_LENGTH / 2)
blend_end    = x_wall[blend_end_i]                                                                # [m]
blend_slope  = np.gradient(trial_thickness, x_wall)[blend_end_i]                                  # [m/m] Trial thickness slope at the end of the blend
blend_len    = blend_end - blend_start                                                            # [m]
s            = (x_wall - blend_start) / blend_len                                                 # 0 to 1 across the blend
blend        = (STRAIGHT_WALL_THICKNESS * (2*s**3 - 3*s**2 + 1) + trial_thickness[blend_end_i] * (3*s**2 - 2*s**3)
                + blend_len * blend_slope * (s**3 - s**2))                                        # [m] Cubic Hermite
thicker_wall = np.where(x_wall <= blend_start, STRAIGHT_WALL_THICKNESS, np.where(x_wall < blend_end, blend, trial_thickness))  # [m]

np.savetxt(os.path.join(output_folder, "trial_wall_thickness.csv"),
           np.column_stack([x_wall, r_inner, iso_depth, trial_thickness, trial_outer_diameter, thicker_wall]), delimiter=",",
           header=f"x [m], inner radius [m], {ISOTHERM_TEMP} K isotherm depth [m], trial wall thickness [m], "
                  f"trial outer diameter [m], thicker wall thickness [m]")

# 1D transient temperature through the thicker wall
depth_frac = np.linspace(0, 1, N_DEPTH)              # fraction of the wall thickness
T_field    = np.zeros((N, N_DEPTH))                  # [K] End of burn temperature field
for i in range(N):
    T_field[i] = wall_transient(h_gas[i], T_aw[i], thicker_wall[i], depth_frac * thicker_wall[i], T_BURN)
T_hot, T_outer, T_mean = T_field[:, 0], T_field[:, -1], np.trapezoid(T_field, depth_frac, axis=1)

#stresses
stress_keys = ["t_th", "t_p", "t_tot", "vm"]
stress_in   = {k: np.zeros(N) for k in stress_keys}   # [Pa] at the hot (inner) face
stress_out  = {k: np.zeros(N) for k in stress_keys}   # [Pa] at the outer face
vm_max      = np.zeros(N)                             # [Pa] max von Mises through the wall
allow_min   = np.zeros(N)                             # [Pa] yield at the worst point / FoS
for i in range(N):
    r = r_inner[i] + depth_frac * thicker_wall[i]
    S = cylinder_stresses(r, T_field[i], wall_pressure[i], youngs_modulus(T_mean[i]), thermal_expansion(T_mean[i]))
    for k in stress_keys:
        stress_in[k][i], stress_out[k][i] = S[k][0], S[k][-1]
    margin       = yield_strength(T_field[i]) / STRESS_FOS / S["vm"]
    worst        = np.argmin(margin)
    vm_max[i]    = S["vm"][worst]
    allow_min[i] = yield_strength(T_field[i][worst]) / STRESS_FOS

yield_hot = yield_strength(T_hot)

# pressure hoop stress, thin vs thick wall
r_over_t    = r_inner / thicker_wall                        # [unitless]
is_thin     = r_over_t >= THIN_WALL_RATIO                   # thin wall valid here
hoop_thin   = wall_pressure * r_inner / thicker_wall        # [Pa] thin wall, p r / t
hoop_thick  = stress_in["t_p"]                              # [Pa] thick wall (Lame), max at the inner face
hoop_design = np.where(is_thin, hoop_thin, hoop_thick)      # [Pa] whichever is valid at each station

# check: series (plane wall) vs Bartz finite volume (cylindrical) solve at the throat with the thicker wall
fv_d, fv_T = bartz.wall_temperature_field(h_gas[throat:throat+1], T_aw[throat:throat+1], r_inner[throat:throat+1],
                                          thicker_wall[throat], STEEL_K, T_BURN)


# plots
figs = {}
IN = c.M2IN

# Trial and thicker wall thickness vs axial position
fig, ax = plt.subplots(); figs["trial_wall_thickness"] = fig
ax.plot(x_rel * IN, thicker_wall * IN, color="black", linewidth=2, label="Thicker wall thickness")
ax.plot(x_rel * IN, iso_depth * IN, color="#d41e11", linewidth=1.5, label=f"{ISOTHERM_TEMP} K isotherm depth")
ax.plot(x_rel * IN, trial_thickness * IN, color="#1f4e9c", linewidth=1.5, linestyle="--", label="Trial wall thickness")
ax.axhline(design_wall * IN, color="gray", linewidth=1, label="Current design wall")
ax.set_xlabel("Axial Position Relative to Throat (in)"); ax.set_ylabel("Thickness (in)")
ax.set_title("Wall Thickness From Isotherm"); ax.grid(True); ax.legend()

# Wall contour (inner + outer surface)
fig, ax = plt.subplots(figsize=(11, 4.5)); figs["trial_wall_contour"] = fig
for s in (1, -1):
    ax.fill_between(x_rel * IN, s * r_inner * IN, s * (r_inner + thicker_wall) * IN, color="gray", alpha=0.5,
                    label="Thicker wall" if s == 1 else None)
    ax.plot(x_rel * IN, s * r_inner * IN, color="red", linestyle="--", label="Inner contour" if s == 1 else None)
    ax.plot(x_rel * IN, s * (r_inner + thicker_wall) * IN, color="black", linewidth=1.5, label="Outer contour" if s == 1 else None)
ax.set_aspect("equal"); ax.grid(True); ax.legend(loc="upper right")
ax.set_xlabel("Axial Position Relative to Throat (in)"); ax.set_ylabel("Radius (in)")
ax.set_title("Chamber Contour With Thicker Wall")

# Temperature field through the thicker wall at end of burn
X_grid = np.repeat(x_rel[:, None], N_DEPTH, axis=1) * IN
R_grid = (r_inner[:, None] + depth_frac[None, :] * thicker_wall[:, None]) * IN
fig, ax = plt.subplots(figsize=(11, 4.5)); figs["trial_wall_temp_field"] = fig
field = ax.pcolormesh(X_grid, R_grid, T_field, shading="gouraud", cmap="inferno", vmin=T_INIT)
ax.contour(X_grid, R_grid, T_field, levels=[ISOTHERM_TEMP], colors="deepskyblue", linewidths=1.5)
ax.plot([], [], color="deepskyblue", label=f"{ISOTHERM_TEMP} K isotherm")
fig.colorbar(field, ax=ax, label="Temperature (K)", shrink=0.8)
ax.set_xlabel("Axial Position Relative to Throat (in)"); ax.set_ylabel("Radius (in)")
ax.set_title(f"1D Transient Wall Temperature, Thicker Wall (t = {T_BURN} s)")
ax.set_aspect("equal"); ax.legend(loc="upper right")

# Hot face / mean / outer face temps
fig, ax = plt.subplots(); figs["trial_wall_surface_temps"] = fig
ax.plot(x_rel * IN, T_hot, color = "#d41e11", label="Hot face")
ax.plot(x_rel * IN, T_mean, color="#e08a00", label="Wall mean")
ax.plot(x_rel * IN, T_outer, color="#1f4e9c", label="Outer face")
ax.axhline(ISOTHERM_TEMP, color="gray", linestyle="--", label=f"{ISOTHERM_TEMP} K")
ax.set_xlabel("Axial Position Relative to Throat (in)"); ax.set_ylabel("Temperature (K)")
ax.set_title("End of Burn Wall Temperatures, Thicker Wall"); ax.grid(True); ax.legend()

# Temperature profile through the wall at the throat over time
times = np.linspace(0, T_BURN, 6)[1:]
d_th  = depth_frac * thicker_wall[throat]
T_t   = wall_transient(h_gas[throat], T_aw[throat], thicker_wall[throat], d_th, times)
fig, ax = plt.subplots(); figs["throat_temp_profile"] = fig
colors = plt.cm.inferno(np.linspace(0.15, 0.8, len(times)))
for T_line, tt, col in zip(T_t, times, colors):
    ax.plot(d_th * IN, T_line, color=col, label=f"t = {tt:.2f} s")
ax.plot(fv_d * IN, fv_T[0], color="black", linestyle=":", label="Finite volume check (cylindrical)")
ax.set_xlabel("Depth Into Wall (in)"); ax.set_ylabel("Temperature (K)")
ax.set_title(f"Throat Wall Temperature Profile (t = {thicker_wall[throat]*IN:.3f} in)"); ax.grid(True); ax.legend()

# Stresses
for key, title, lines in [
    ("thermal_stress",  "Thermal Hoop Stress, Thicker Wall",
        [(stress_in["t_th"], "Hot face", "#d41e11"), (stress_out["t_th"], "Outer face", "#1f4e9c")]),
    ("pressure_stress", "Pressure Hoop Stress, Thicker Wall",
        [(stress_in["t_p"], "Hot face", "#d41e11"), (stress_out["t_p"], "Outer face", "#1f4e9c")]),
    ("total_hoop_stress", "Total Hoop Stress, Thicker Wall",
        [(stress_in["t_tot"], "Hot face", "#d41e11"), (stress_out["t_tot"], "Outer face", "#1f4e9c")]),
]:
    fig, ax = plt.subplots(); figs[key] = fig
    for y, lab, col in lines:
        ax.plot(x_rel * IN, y * c.PA2PSI, color=col, linewidth=1, label=lab)
    ax.plot(x_rel * IN, yield_hot * c.PA2PSI, color=(0, 0.3, 0.7), linestyle="--", linewidth=1, label="Yield")
    ax.plot(x_rel * IN, -yield_hot * c.PA2PSI, color=(0, 0.3, 0.7), linestyle="--", linewidth=1)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_xlabel("Axial Position Relative to Throat (in)"); ax.set_ylabel("Stress (psi)")
    ax.set_title(title); ax.grid(True); ax.legend(loc="lower left")

# Pressure hoop stress, thin vs thick wall
fig, ax = plt.subplots(); figs["hoop_stress_thin_vs_thick"] = fig
ax.fill_between(x_rel * IN, 0, 1, where=~is_thin, transform=ax.get_xaxis_transform(), color="gray", alpha=0.2,
                label=f"Thick wall (r/t < {THIN_WALL_RATIO})")
ax.plot(x_rel * IN, hoop_thin * c.PA2PSI, color="#e08a00", linewidth=1, linestyle="--", label="Thin wall (p r / t)")
ax.plot(x_rel * IN, hoop_thick * c.PA2PSI, color="#1f4e9c", linewidth=1, linestyle="--", label="Thick wall (Lame, inner face)")
ax.plot(x_rel * IN, hoop_design * c.PA2PSI, color="black", linewidth=1.5, label="Valid model")
ax.set_xlabel("Axial Position Relative to Throat (in)"); ax.set_ylabel("Stress (psi)")
ax.set_title("Pressure Hoop Stress, Thicker Wall"); ax.grid(True); ax.legend(loc="upper left")

ax2 = ax.twinx()
ax2.plot(x_rel * IN, r_over_t, color="gray", linewidth=0.8, linestyle=":")
ax2.axhline(THIN_WALL_RATIO, color="gray", linewidth=0.5)
ax2.set_ylabel("r / t (dotted)")


# printing stuff
print(f"Design wall: {design_wall*IN:.3f} in")
print(f"Max {ISOTHERM_TEMP} K isotherm depth: {iso_depth.max()*IN:.3f} in (throat: {iso_depth[throat]*IN:.3f} in)")
print(f"Trial thickness range: {trial_thickness.min()*IN:.3f} to {trial_thickness.max()*IN:.3f} in "
      f"(floor {MIN_THICKNESS*IN:.3f} in at {np.sum(iso_depth*THICKNESS_FOS < MIN_THICKNESS)} stations)")
if capped.any():
    print(f"WARNING: isotherm reached the full design wall at {capped.sum()} stations, trial thickness is capped there")
print(f"\nThicker wall, end of burn: max hot face {T_hot.max():.0f} K, max mean {T_mean.max():.0f} K, max outer face {T_outer.max():.0f} K")
print(f"Throat check, hot face: series (infinite cylinder) {T_field[throat, 0]:.0f} K vs finite volume (annular wall) {fv_T[0, 0]:.0f} K")
print(f"\nMax thermal hoop stress magnitude: {np.abs(stress_in['t_th']).max()*c.PA2PSI:.0f} psi (hot face, compressive)")
print(f"Thin wall valid (r/t >= {THIN_WALL_RATIO}) at {is_thin.sum()} of {N} stations, thick wall at {(~is_thin).sum()}"
      f" (throat r/t = {r_over_t[throat]:.1f}, {'thin' if is_thin[throat] else 'thick'})")
if (~is_thin).any():
    print(f"  thick wall region: {x_rel[~is_thin].min()*IN:.2f} to {x_rel[~is_thin].max()*IN:.2f} in from throat")
print(f"Max pressure hoop stress: thin wall {hoop_thin.max()*c.PA2PSI:.0f} psi, thick wall {hoop_thick.max()*c.PA2PSI:.0f} psi, "
      f"valid model {hoop_design.max()*c.PA2PSI:.0f} psi")
print(f"Throat pressure hoop stress: thin wall {hoop_thin[throat]*c.PA2PSI:.0f} psi, thick wall {hoop_thick[throat]*c.PA2PSI:.0f} psi")
print(f"Min von Mises margin (target / von Mises): {np.min(allow_min / vm_max):.2f}")
fails = vm_max > allow_min
if fails.any():
    print(f"WARNING: von Mises exceeds yield/{STRESS_FOS} at {fails.sum()} stations "
          f"({x_rel[fails].min()*IN:.2f} to {x_rel[fails].max()*IN:.2f} in from throat)")


#export plots
for name, fig in figs.items():
    fig.tight_layout()
    fig.savefig(os.path.join(output_folder, f"{name}.png"), dpi=200)

#plt.show()
