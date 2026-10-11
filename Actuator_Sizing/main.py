import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
















import matplotlib.pyplot as plt
plt.style.use('dark_background')
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
# OUTPUTS on or off
show_outputs = False
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

Cv_regulator = 0.06 # Cv for BCLS N2 Regulator (Should be the bottleneck)
dead_volume_m3 = 0.5 * c.IN2M**3  # PLACEHOLDER - replace with actual tubing+fitting+clearance volume in in^3
T_ambient_R = 530  # PLACEHOLDER - 70F, replace if you know actual ambient/supply gas temp

Cv_accumulator = 2  # PLACEHOLDER - Cv of accumulator outlret
accumulator_volume_m3 = 25 * c.IN2M**3  # PLACEHOLDER - accumulator volume in in^3

Cv_solenoid = 1

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

def sum_Cv(*Cvs):
    total_Cv = 1 / np.sqrt(sum(1 / cv**2 for cv in Cvs))
    if show_outputs == True:
        print(f"Total Cv: {total_Cv:.2f}")
    return total_Cv

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

def actuation_time_kinematics_accumulator(rod_mass, piston_diameter, arm_length, friction_total, force_valve, total_Cv, fill_pressure_psig, accumulator_volume_m3, dead_volume_m3, show_outputs, T_ambient_R=530, gas_SG=0.967, source_Cv=0.0):
    R_SPECIFIC_N2 = 296.8  # J/(kg*K)
    P_ATM_PSIA = 14.7
    P_ATM_PA = P_ATM_PSIA * c.PSI2PA
    T_STD_K = 288.71  # 60 F, standard reference temp for SCFH
    density_std = (P_ATM_PSIA * c.PSI2PA) / (R_SPECIFIC_N2 * T_STD_K)  # kg/m^3 of N2 at standard conditions

    P_fill_psia = fill_pressure_psig + P_ATM_PSIA
    T_K = T_ambient_R * 5 / 9  # Rankine -> Kelvin
    piston_area = np.pi * (piston_diameter / 2)**2

    # Chamber starts vented to atmosphere, accumulator starts fully charged
    chamber_mass = (P_ATM_PA * dead_volume_m3) / (R_SPECIFIC_N2 * T_K)
    accumulator_mass = (P_fill_psia * c.PSI2PA * accumulator_volume_m3) / (R_SPECIFIC_N2 * T_K)
    accumulator_mass_full = accumulator_mass  # fully charged mass; a connected source can refill back up to this but no further
    source_mass_total = 0.0  # running total of gas the regulator has put back into the accumulator
    chamber_mass_gained = 0.0  # running total of gas moved accumulator -> chamber

    time = 0
    time_step = 0.00001
    max_time = 5.0  # cutoff in case the force balance never reaches 90 deg
    piston_velocity = 0
    dist_travelled = 0
    valve_angle = 0
    stalled = False
    time_history = []
    angle_history = []
    velocity_history = []
    volume_swept_history = []
    distance_travelled_history = []
    chamber_pressure_history = []
    accumulator_pressure_history = []
    flow_scfm_history = []
    source_flow_scfm_history = []

    while valve_angle <= 90 and time < max_time:
        V_chamber = dead_volume_m3 + piston_area * dist_travelled
        P_chamber_Pa = chamber_mass * R_SPECIFIC_N2 * T_K / V_chamber
        P_chamber_psia = P_chamber_Pa / c.PSI2PA
        P_acc_psia = accumulator_mass * R_SPECIFIC_N2 * T_K / accumulator_volume_m3 / c.PSI2PA

        # Upstream pressure is now the (falling) accumulator pressure instead of a fixed supply
        if P_chamber_psia < P_acc_psia:
            if P_acc_psia >= 2 * P_chamber_psia:
                Q_scfh = 816 * total_Cv * P_acc_psia / np.sqrt(gas_SG * T_ambient_R)  # choked
            else:
                Q_scfh = 962 * total_Cv * np.sqrt((P_acc_psia**2 - P_chamber_psia**2) / (gas_SG * T_ambient_R))  # subsonic
        else:
            Q_scfh = 0
        Q_scfm = Q_scfh / 60
        mdot = (Q_scfh / 3600) * 0.0283168 * density_std  # kg/s (0.0283168 = m^3 per ft^3)

        # Mass leaves the accumulator and enters the chamber. Cap the transfer at the amount that would just equalize
        # the two pressures, so a big Cv with a tiny dead volume can't overshoot in a single time step.
        mass_chamber_equalized = (accumulator_mass + chamber_mass) * V_chamber / (V_chamber + accumulator_volume_m3)
        dm = min(mdot * time_step, max(mass_chamber_equalized - chamber_mass, 0))
        chamber_mass += dm
        accumulator_mass -= dm
        chamber_mass_gained += dm

        # Regulator/source still connected to the accumulator: it refills the accumulator through source_Cv.
        # Upstream is the regulator outlet (held at P_fill_psia), downstream is the accumulator, so this is the same
        # choked/subsonic Cv equation with the accumulator as the 'chamber'. Set source_Cv=0 for an isolated accumulator.
        Q_src_scfh = 0
        if source_Cv > 0 and P_acc_psia < P_fill_psia:
            if P_fill_psia >= 2 * P_acc_psia:
                Q_src_scfh = 816 * source_Cv * P_fill_psia / np.sqrt(gas_SG * T_ambient_R)  # choked
            else:
                Q_src_scfh = 962 * source_Cv * np.sqrt((P_fill_psia**2 - P_acc_psia**2) / (gas_SG * T_ambient_R))  # subsonic
            mdot_src = (Q_src_scfh / 3600) * 0.0283168 * density_std  # kg/s
            dm_src = min(mdot_src * time_step, max(accumulator_mass_full - accumulator_mass, 0))  # can't overfill past the original charge
            accumulator_mass += dm_src
            source_mass_total += dm_src

        chamber_force = max(P_chamber_Pa - P_ATM_PA, 0) * piston_area
        f_net = chamber_force - force_valve - friction_total
        if piston_velocity <= 0 and f_net <= 0:
            f_net = 0  # not enough pressure built up yet to overcome friction/valve resistance
            # Stalled = chamber has equalized with the accumulator, and (if a source is connected) the accumulator has
            # also recharged all the way to regulator pressure, so there is nothing left to build more pressure, and the piston still can't move.
            if P_chamber_psia >= P_acc_psia - 0.01 and (source_Cv == 0 or P_acc_psia >= P_fill_psia - 0.01):
                stalled = True
                break

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
        accumulator_pressure_history.append(P_acc_psia)
        flow_scfm_history.append(Q_scfm)
        source_flow_scfm_history.append(Q_src_scfh / 60)
        time += time_step

    reached_90 = valve_angle > 90
    if stalled:
        print(f"WARNING: piston stalled at t={time:.3f}s, accumulator and chamber equalized at {P_acc_psia - P_ATM_PSIA:.1f} psig — accumulator volume too small (or total_Cv/force balance issue)")
    elif not reached_90:
        print(f"WARNING: valve did not reach 90 deg within {max_time}s cutoff — check accumulator volume/Cv/force balance")
    if not reached_90:
        time = np.inf

    if show_outputs == True:
        plt.figure(figsize=(14, 8))
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
        plt.ylim(bottom = min(np.min(velocity_history) - (0.05 * (np.max(velocity_history) - np.min(velocity_history))), 0))
        plt.grid(True)

        plt.subplot(2, 3, 3)
        plt.plot(time_history, distance_travelled_history)
        plt.xlabel("Time [s]")
        plt.ylabel("Distance Swept [m]")
        plt.title("Time vs Distance Swept")
        plt.grid(True)

        plt.subplot(2, 3, 4)
        plt.plot(time_history, chamber_pressure_history)
        plt.axhline(P_fill_psia, color='r', linestyle='--', label='Accumulator fill')
        plt.ylim(bottom=0)
        plt.xlabel("Time [s]")
        plt.ylabel("Chamber Pressure [psia]")
        plt.title("Chamber Pressure Fill")
        plt.legend()
        plt.grid(True)

        plt.subplot(2, 3, 5)
        plt.plot(time_history, flow_scfm_history, label='Accumulator discharge')
        if source_Cv > 0:
            plt.plot(time_history, source_flow_scfm_history, linestyle='--', label='Regulator refill')
            plt.legend()
        plt.ylim(bottom=0)
        plt.xlabel("Time [s]")
        plt.ylabel("Flow [SCFM]")
        plt.title("Accumulator Discharge Flow")
        plt.grid(True)

        plt.subplot(2, 3, 6)
        plt.plot(time_history, accumulator_pressure_history, label='Accumulator')
        plt.plot(time_history, chamber_pressure_history, label='Actuator chamber', linestyle=':')
        plt.ylim(bottom=0)
        plt.xlabel("Time [s]")
        plt.ylabel("Pressure [psia]")
        plt.title("Accumulator Blowdown")
        plt.legend()
        plt.grid(True)

        plt.tight_layout()
        plt.show()

        if reached_90:
            print(f"Actuation time (accumulator-fed): {time:.3f}s")
        else:
            print("Actuation time (accumulator-fed): valve never reached 90 deg")
        print(f"Peak accumulator discharge flow: {max(flow_scfm_history):.2f} SCFM")
        if source_Cv > 0:
            print(f"Peak regulator refill flow: {max(source_flow_scfm_history):.2f} SCFM, {100 * source_mass_total / max(chamber_mass_gained, 1e-12):.1f}% as much gas as the accumulator delivered to the chamber")
        print(f"Final accumulator pressure: {accumulator_pressure_history[-1]:.1f} psia ({accumulator_pressure_history[-1] - P_ATM_PSIA:.1f} psig) from {fill_pressure_psig:.0f} psig fill")
        print(f"Final chamber pressure: {chamber_pressure_history[-1]:.1f} psia ({chamber_pressure_history[-1] - P_ATM_PSIA:.1f} psig)")

    return volume_swept_history, time_history, angle_history, time, chamber_pressure_history, flow_scfm_history, accumulator_pressure_history


piston_force = pressure * np.pi * ((piston_diameter**2) / 4)
if show_outputs == True:
    print(f'Maximum possible net force disregarding friction (and valve arm if real condition): {piston_force * c.N2LBF:.2f}')

required_torque, arm_length, torque = calc_torque_piston(breaking_torque, breaking_torque_safety_factor, piston_force, piston_stroke_length, show_outputs)
force_valve = breaking_torque * np.sqrt(2) / arm_length
friction_total = calc_friction_force(piston_seal_length, shaft_seal_length, piston_seal_area, shaft_seal_area, 0)

Cv_total = sum_Cv(Cv_accumulator, Cv_solenoid)

#volume_swept_history, time_history, angle_history, time, chamber_pressure_history, flow_scfm_history = actuation_time_kinematics_flow_limited(
#        rod_mass, piston_diameter, arm_length, friction_total, force_valve, Cv, pressure * c.PA2PSI, dead_volume_m3, show_outputs, T_ambient_R)
if __name__ == "__main__":
    acc_volume_swept_history, acc_time_history, acc_angle_history, acc_time, acc_chamber_pressure_history, acc_flow_scfm_history, accumulator_pressure_history = actuation_time_kinematics_accumulator(
        rod_mass, piston_diameter, arm_length, friction_total, force_valve, Cv_total, pressure * c.PA2PSI, accumulator_volume_m3, dead_volume_m3, show_outputs, T_ambient_R, source_Cv=Cv_regulator)