import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
















import matplotlib.pyplot as plt
import numpy as np

# sorry
import sys,pathlib,collections,importlib.abc; r=next(p for p in pathlib.Path(__file__).resolve().parents if p.name=="PSPL_Rocket_A");   m=collections.defaultdict(list); [m[f.stem.casefold()].append(f) for f in r.rglob("*.py") if f.stem.isidentifier() and f.name != "__init__.py" and not (set(f.relative_to(r).parts)&{".git",".venv","__pycache__","build","dist"})]; dup={k:v for k,v in m.items() if   len(v)>1}; sys.meta_path.insert(0,type("AmbiguousBareImportBlocker",(importlib.abc.MetaPathFinder,),{"find_spec":lambda   self,fullname,path=None,target=None: (_ for _ in ()).throw(ImportError(f"The import {fullname!r} could refer to any following packages: "+" ".join("\n\t" + str(p.relative_to(r)) for p in dup[fullname.casefold()])+f"\n\nSpecify which package it is by using the folder.\nFor example:\n\t'import SFD.{fullname}'")) if "." not in fullname and fullname.casefold() in dup else None})()); sys.path.insert(0,str(r)) if str(r) not   in sys.path else None; [sys.path.append(str(v[0].parent)) for k,v in m.items() if k not in dup and str(v[0].parent) not in sys.path]

import constants as c # type: ignore
import vehicle_parameters_functions # type: ignore
import vehicle_parameters # type: ignore
import vehicle_main # type: ignore
import print_filter # type: ignore

# def Name_of_Script_Main_Function(parameters): # this is a separate function for reasons i dont remember but its important...
    





















###############################
# OUTPUTS "on" (1) or "off" (0)
show_outputs = True
###############################

###############################
# INPUT PARAMETERS
###############################
breaking_torque = 242 * c.IN_LB2NM
breaking_torque_safety_factor = 3
piston_stroke_length = 4 * c.IN2M
rod_mass = 1.5456 * c.LBM2KG # Estimated from CAD
piston_diameter = 2.5 * c.IN2M
piston_retracted_length = 9.44 * c.IN2M
piston_extended_length = piston_retracted_length + piston_stroke_length
shaft_diameter = 0.625 * c.IN2M
piston_seal_length = np.pi * piston_diameter
shaft_seal_length = np.pi * shaft_diameter
piston_seal_area = 0.21 * c.IN2M * piston_seal_length # worst case scenario, 300 series
shaft_seal_area = 0.21 * c.IN2M * shaft_seal_length
pressure = 100 * c.PSI2PA

def calc_friction_force(piston_seal_length, shaft_seal_length, piston_seal_area, shaft_seal_area, show_outputs):
    fc_piston = (4 * piston_seal_length * c.M2IN) * c.LBF2N # assuming 4, worst case, for now
    fc_shaft = (4 * shaft_seal_length * c.M2IN) * c.LBF2N
    fh_piston = (18 * piston_seal_area * c.M22IN2) * c.LBF2N # from parker oring handbook figure 5-10, lowkey using it for u-cup :skull:
    fh_shaft = (18 * shaft_seal_area * c.M22IN2) * c.LBF2N
    friction_piston = fc_piston + fh_piston
    friction_shaft = fc_shaft + fh_shaft
    total_friction = friction_piston * 2 + friction_shaft # two seals on piston, 1 on rod
    if show_outputs == True:
        print(f"fc_piston: {fc_piston * c.N2LBF:.2f} LBF")
        print(f"fc_shaft: {fc_shaft * c.N2LBF:.2f} LBF")
        print(f"fh_piston: {fh_piston * c.N2LBF:.2f} LBF")
        print(f"fh_shaft: {fh_shaft * c.N2LBF:.2f} LBF")
        print(f"F_total_friction: {total_friction * c.N2LBF:.2f} LBF")
    return total_friction
    
def calc_torque_piston(breaking_torque, breaking_torque_safety_factor, piston_force, piston_stroke_length, show_outputs):
    required_torque = breaking_torque * breaking_torque_safety_factor
    arm_length = piston_stroke_length / np.sqrt(2)
    torque = arm_length * piston_force / np.sqrt(2)
    if show_outputs == True:
        print(f"The piston will produce ~{torque * c.NM2IN_LB:.2f} lb-in torque at {pressure * c.PA2PSI} psi.")
        print(f"The required torque with a safety factor of 3 is {required_torque * c.NM2IN_LB:.2f}")
        print(f"Length of valve arm would be {arm_length * c.M2IN:.2f}")
    return required_torque, arm_length, torque

def actuation_time_kinematics_flow_limited(rod_mass, piston_diameter, arm_length, friction_total, force_valve, Cv, supply_pressure_psig, dead_volume_m3, show_outputs, T_ambient_R=530, gas_SG=0.967):
    R_SPECIFIC_N2 = 296.8  # J/(kg*K)
    P_ATM_PSIA = 14.7
    P_ATM_PA = P_ATM_PSIA * c.PSI2PA
    T_STD_K = 288.71  # 60 F, standard reference temp for SCFH
    density_std = (P_ATM_PSIA * c.PSI2PA) / (R_SPECIFIC_N2 * T_STD_K)  # kg/m^3 of N2 at standard conditions

    P_supply_psia = supply_pressure_psig + P_ATM_PSIA
    T_K = T_ambient_R * 5 / 9  # Rankine -> Kelvin (same absolute zero, just rescaled)
    piston_area = np.pi * (piston_diameter / 2)**2

    # Chamber starts vented to atmosphere, not vacuum
    chamber_mass = (P_ATM_PA * dead_volume_m3) / (R_SPECIFIC_N2 * T_K)
 
    time = 0
    time_step = 0.00001
    max_time = 5.0  # cutoff in case Cv/force balance never reaches 90 deg
    piston_velocity = 0
    dist_travelled = 0
    valve_angle = 0
    time_history = []
    angle_history = []
    velocity_history = []
    volume_swept_history = []
    distance_travelled_history = []
    chamber_pressure_history = []
    flow_scfm_history = []
 
    while valve_angle <= 90 and time < max_time:
        V_chamber = dead_volume_m3 + piston_area * dist_travelled
        P_chamber_Pa = chamber_mass * R_SPECIFIC_N2 * T_K / V_chamber
        P_chamber_psia = P_chamber_Pa / c.PSI2PA
 
        if P_chamber_psia < P_supply_psia:
            if P_supply_psia >= 2 * P_chamber_psia:
                Q_scfh = 816 * Cv * P_supply_psia / np.sqrt(gas_SG * T_ambient_R)  # choked
            else:
                Q_scfh = 962 * Cv * np.sqrt((P_supply_psia**2 - P_chamber_psia**2) / (gas_SG * T_ambient_R))  # subsonic
        else:
            Q_scfh = 0
        Q_scfm = Q_scfh / 60
        mdot = (Q_scfh / 3600) * 0.0283168 * density_std  # kg/s (0.0283168 = m^3 per ft^3)
        chamber_mass += mdot * time_step
 
        chamber_force = max(P_chamber_Pa - P_ATM_PA, 0) * piston_area
        f_net = chamber_force - force_valve - friction_total
        if piston_velocity <= 0 and f_net <= 0:
            f_net = 0  # not enough pressure built up yet to overcome friction/valve resistance
 
        dist_travelled = max(dist_travelled + piston_velocity * time_step + 0.5 * (f_net / rod_mass) * time_step**2, 0)
        piston_velocity = max(piston_velocity + (f_net / rod_mass) * time_step, 0)
 
        if dist_travelled == 0:
            valve_angle = 0
        else:
            valve_angle = np.degrees((np.pi / 2) - np.arctan((((arm_length / dist_travelled) - (1 / np.sqrt(2))) * np.sqrt(2))))
        volume_swept = dist_travelled * piston_area
 
        time_history.append(time)
        angle_history.append(valve_angle)
        velocity_history.append(piston_velocity)
        volume_swept_history.append(volume_swept)
        distance_travelled_history.append(dist_travelled)
        chamber_pressure_history.append(P_chamber_psia)
        flow_scfm_history.append(Q_scfm)
        time += time_step
 
    if time >= max_time:
        print(f"WARNING: valve did not reach 90 deg within {max_time}s cutoff — check Cv/force balance")
 
    if show_outputs == True:
        plt.subplot(2, 3, 1)
        plt.plot(time_history, angle_history)
        plt.xlabel("Time [s]")
        plt.ylabel("Valve Angle [˚]")
        plt.title("Valve Angle Over Time")
        plt.ylim(bottom = -5, top = 95)
        plt.yticks(np.arange(0, plt.ylim()[1] + 10, 10))
        plt.grid(True)
 
        plt.subplot(2, 3, 2)
        plt.plot(time_history, velocity_history)
        plt.xlabel("Time [s]")
        plt.ylabel("Velocity [m/s]")
        plt.title("Time vs Velocity")
        plt.ylim(bottom = np.min(np.min(velocity_history) - (0.05 * (np.max(velocity_history) - np.min(velocity_history))), 0))
        plt.grid(True)
 
        plt.subplot(2, 3, 3)
        plt.plot(time_history, distance_travelled_history)
        plt.xlabel("Time [s]")
        plt.ylabel("Distance Swept [m]")
        plt.title("Time vs Distance Swept")
        plt.grid(True)
 
        plt.subplot(2, 3, 4)
        plt.plot(time_history, chamber_pressure_history)
        plt.axhline(P_supply_psia, color='r', linestyle='--', label='Supply')
        plt.ylim(bottom=0)
        plt.xlabel("Time [s]")
        plt.ylabel("Chamber Pressure [psia]")
        plt.title("Chamber Pressure Fill")
        plt.legend()
        plt.grid(True)
 
        plt.subplot(2, 3, 5)
        plt.plot(time_history, flow_scfm_history)
        plt.ylim(bottom=0)
        plt.xlabel("Time [s]")
        plt.ylabel("Flow [SCFM]")
        plt.title("Solenoid Flow Demand")
        plt.grid(True)

        plt.tight_layout()
        plt.show()
 
        print(f"Actuation time (flow-limited): {time:.3f}s")
        print(f"Peak solenoid flow demand: {max(flow_scfm_history):.2f} SCFM")
        print(f"Final chamber pressure: {chamber_pressure_history[-1]:.1f} psia ({chamber_pressure_history[-1] - P_ATM_PSIA:.1f} psig) vs supply {supply_pressure_psig:.0f} psig")
 
    return volume_swept_history, time_history, angle_history, time, chamber_pressure_history, flow_scfm_history

piston_force = pressure * np.pi * ((piston_diameter**2) / 4)
if show_outputs == True:
    print(f'Maximum possible net force disregarding friction (and valve arm if real condition): {piston_force * c.N2LBF:.2f}')

required_torque, arm_length, torque = calc_torque_piston(breaking_torque, breaking_torque_safety_factor, piston_force, piston_stroke_length, show_outputs)
force_valve = breaking_torque * np.sqrt(2) / arm_length
friction_total = calc_friction_force(piston_seal_length, shaft_seal_length, piston_seal_area, shaft_seal_area, 0)

# Constants
Cv = 0.06 # Cv for BCLS N2 Regulator (Should be the bottleneck)
dead_volume_m3 = 0.5 * c.IN2M**3  # PLACEHOLDER - replace with actual tubing+fitting+clearance volume in in^3
T_ambient_R = 530  # PLACEHOLDER - 70F, replace if you know actual ambient/supply gas temp

volume_swept_history, time_history, angle_history, time, chamber_pressure_history, flow_scfm_history = actuation_time_kinematics_flow_limited(
        rod_mass, piston_diameter, arm_length, friction_total, force_valve, Cv, pressure * c.PA2PSI, dead_volume_m3, show_outputs, T_ambient_R)