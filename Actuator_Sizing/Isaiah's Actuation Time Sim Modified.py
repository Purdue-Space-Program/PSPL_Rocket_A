#   Written by Isaiah Jarvis
#
#   This code runs a simulation of actuation time
#
#   Variables are iterated over and plotted to show how they could affect our design
# 
#   If you want to find the actuation time of a single model simply delete/comment out the
#   graphing section (lines 56 - 184), update the constants, and run the "actation_time" function
# 


from fluids import core
from fluids import fittings
import math as m
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
from pint import UnitRegistry
import time
u = UnitRegistry()

# Tank Pressure / fuel properties
tank_pressure = 287.0 * u.psi                    # Ethanol tank pressure that will provide power to the actuator
ambient_pressure = 14.696 * u.psi                # Fluid in the back of the piston is thrown off-board. Ambient pressure acts as that "Tank pressure"
temp = 173.0 * u.kelvin                          # Approximate ethanol temperature in the tank 293 is ambient

# Calculated properties
# Fluid properties - Used to calculate pressure drops
density = CP.PropsSI('D', 'T', temp.magnitude, 'P', tank_pressure.to('Pa').magnitude, 'ETHANOL') * (u.kg / u.m**3)
dynamic_mu = CP.PropsSI('V', 'P', tank_pressure.magnitude, 'T', temp.magnitude, 'ETHANOL') * (u.pascal * u.second) #viscosity
kinetic_mu = dynamic_mu / density

print(kinetic_mu.to(u.cSt))


# Start of analysis
def main():
    # Tank Pressure / fuel properties
    tank_pressure = 250.0 * u.psi                    # Has to be defined in the main function to be used in loops

    #Tube constants
    tube_outer_diameter = 0.25 * u.inch
    tube_thickness = 0.035 * u.inch

    # Piston geometry
    piston_radius = 1 * u.inch                     
    piston_mass = (100 * u.gram).to(u.kg)            # This value can be imported from CAD - Estimated for now

    # Valve constants
    valve_arm = 1.4 * u.inch                         # Length of valve arm
    valve_torque = 120 * u.inch * u.lbf              # Force required to turn the valve - Using CMS values
    actuation_angle = (90 * u.degree).to(u.radian)   # Angle required to open/close the valve

    # Shaft geometry
    shaft_length = 4 * u.inch                        # Length from actuator hinge to shaft-valve arm hinge
    shaft_diameter = 0.32 * u.inch
    theta = (135 * u.degree).to(u.radian)            # Angle between shaft arm and valve arm

    # O-ring constants
    rotary_seal_length = 0.78                        # Length of O-ring that makes contact with seal surface
    rotary_seal_area = 0.05                          # Area of O-ring that makes contact with seal surface

    piston_seal_length = 3.14
    piston_seal_area = 0.17

    # Solenoid Constants
    # Cv = 3.0                                         # Will come from solenoid manufacturer
    # solenoid_diameter = tube_outer_diameter          # Most likely the same as tube size - Will also come from manufacturer

    
    test = manifold_dp(35 * u.inch**3 / u.second)
    test2 = solenoid_dp(35 * u.inch**3 / u.second, 32)
    print(test)
    print(test2)

    forward_time = actuation_time(piston_mass, tube_outer_diameter, tube_thickness, piston_radius, valve_arm, valve_torque, shaft_length, shaft_diameter, theta, actuation_angle, tank_pressure, rotary_seal_length, rotary_seal_area, piston_seal_length, piston_seal_area)

    print(forward_time)

def actuation_time(piston_mass, tube_outer_diameter, tube_thickness, piston_radius, valve_arm, valve_torque, shaft_length, shaft_diameter, theta, actuation_angle, tank_pressure, rotary_seal_length, rotary_seal_area, piston_seal_length, piston_seal_area):
    # Geometry calculations
    piston_area = m.pi * piston_radius ** 2
    shaft_area = m.pi * (shaft_diameter / 2) ** 2

    #Triangle properties
    hinge_distance = m.sqrt((shaft_length**2 + valve_arm**2 - 2 * shaft_length * valve_arm * m.cos(theta.magnitude)).magnitude) * u.inch    # Distance between the actuator hinge and valve - static value
    valve_angle_start = m.acos((-(shaft_length**2) + valve_arm**2 + hinge_distance**2) / ((2) * (valve_arm) * (hinge_distance))) * u.radian       # Angle between hinge points and valve arm
    valve_angle = valve_angle_start

    # Tubes
    tube_area = m.pi * (tube_outer_diameter / 2 - tube_thickness)**2

    #force of friction from o-ring from parker guide found here: https://www.parker.com/literature/O-Ring%20Division%20Literature/ORD%205700%20Parker_O-Ring_Handbook.pdf
    # OLD FRICTION CALCULATION
    # lp = m.pi #length of the o ring (in)
    # fc = 1 #Friction due to compression of o ring (lbf)
    # fh = 25 #total friction due to hydraulic pressure on seal (lb)
    # a_p = (1/16) * m.pi #projected area of the o ring (in)
    # FC = lp * fc #Friction due to seal compression of the o ring (lbf)
    # FH = fh * a_p #Hydraulic pressure on o ring around piston (lbf)
    # F_friction = ((FC + FH) * u.lbf).to(u.newton) #force of friction on the o ring around piston (lbf)

    # Friction due to hydraulic pressure will be calculated in the loop as the pressure changes
    Fc_rotary = ((2.4 * rotary_seal_length) * u.lbf).to(u.newton) # Friction due to the sealing surface on rotarty seals
    Fc_piston = ((2.4 * piston_seal_length) * u.lbf).to(u.newton) # Friction due to the sealing surface on piston seal
    

    forward_t = 0 * u.second                        # t = 0 is when the solenoid is opened
    time_step = 0.00004 * u.second                   # Each step of time

    piston_velocity = 0 * u.inch / u.second         # Will be updated
    distance_traveled = 0 * u.inch                  # Will be updated

    while valve_angle.magnitude <= (valve_angle_start + actuation_angle).magnitude:
        forward_volumetric_flow = piston_velocity * piston_area
        backward_volumetric_flow = piston_velocity * (piston_area - shaft_area)

        forward_line_velocity = forward_volumetric_flow / tube_area     # As piston speeds up the velocity required to keep pushing will increase
        backward_line_velocity = backward_volumetric_flow / tube_area

        #These variables will be used to calculate the net force on the piston
        pressure_full = (tank_pressure - 2 * bend_dp(forward_line_velocity, 90 * u.degree, 0.5 * u.inch, tube_outer_diameter, tube_thickness, 0.0015 * u.mm) - solenoid_dp(forward_volumetric_flow, 34) - manifold_dp(forward_volumetric_flow))
        pressure_shaft = ambient_pressure + 2 * bend_dp(backward_line_velocity, 90 * u.degree, 0.5 * u.inch, tube_outer_diameter, tube_thickness, 0.0015 * u.mm) + solenoid_dp(backward_volumetric_flow, 21) + manifold_dp(backward_volumetric_flow)
        
        F_full = pressure_full * piston_area
        F_shaft = pressure_shaft * (piston_area - shaft_area)
        F_valve_resistance = valve_torque / (valve_arm * m.sin(theta.magnitude))

        Fh_piston= ((10 + 0.03 * abs((pressure_full - pressure_shaft).magnitude)) * piston_seal_area * u.lbf).to(u.newton)
        Fh_rotary = ((10 + 0.03 * pressure_shaft.magnitude) * rotary_seal_area * u.lbf).to(u.newton)
        #Fh_piston = (fh * piston_seal_area * u.lbf).to(u.newton)
        F_rotary = Fc_rotary + Fh_rotary
        F_piston = Fc_piston + Fh_piston

        F_net = (F_full - F_shaft - 2*F_rotary - F_piston - F_valve_resistance).to(u.newton)

        # Possibly necessary print statements
        #print(f"Time = {forward_t:.4f}")
        # print(f"Theta: {theta.to(u.degree):.2f}")
        # print(f"Start Valve Angle: {valve_angle_start.to(u.degree):.2f}")
        # print(f"Valve Angle: {valve_angle.to(u.degree):.2f}")
        # print(f"line velocity: {forward_line_velocity:.4f}")
        # print(f"Upstream force: {F_full:2f}")
        #print(f"Upstream pressure: {pressure_full:.2f}")
        # print(f"Downstream force: {F_shaft:2f}")
        #print(f"Downstream pressure: {pressure_shaft:.2f}")
        # print(f"Friction force: {2 * F_rotary + F_piston:2f}")
        # print(f"Valve resistance force: {F_valve_resistance:2f}")
        # print(f"Net Force: {F_net:.2f}")
        # print(f"Distance traveled: {distance_traveled:.4f}")
        # print(f"Volumetric Flow: {forward_volumetric_flow}")
        #print(f"Piston Velocity 1: {piston_velocity:.5f}")
        # print(F_rotary)
        # print(F_piston)
        # time.sleep(0.25)
        #print()

        piston_velocity += (F_net * time_step) / (piston_mass)
        distance_traveled += piston_velocity * time_step + ((F_net * (time_step ** 2)) / (2 * piston_mass))

        # Law of Cosine
        theta = m.acos(
            (((distance_traveled + shaft_length)**2 + valve_arm**2 - hinge_distance**2) / (2 * (distance_traveled + shaft_length) * valve_arm)).magnitude
        ) * u.radian

        valve_angle = valve_angle = m.acos((-((distance_traveled + shaft_length)**2) + valve_arm**2 + hinge_distance**2) / ((2) * (valve_arm) * (hinge_distance))) * u.radian

        forward_t += time_step

    piston_velocity = 0 * u.inch / u.second     # reset the velocity for back actuation

    return forward_t

# These equations are from manufacturer specifications using hydraulic oil. In turbulent flow, the viscosities of the fluids don't differ that much.
# Parker DSH104 Solenoid Spool
def solenoid_dp(vol_flow, ports):
    vol_flow = vol_flow.to(u.gallon / u.minute).magnitude

    # DSH104
    # if ports == 32 or 34:
    #     dp = (1.7142 * vol_flow**2 + 8.7857 * vol_flow) * u.psi
    # else:
    #     dp = (2.7143 * vol_flow**2 - 2.1143 * vol_flow) * u.psi

    #DSH084
    dp = (3.2857 * vol_flow**2 + 4.8571 * vol_flow) * u.psi

    return dp


# Parker B10-4-8T Manifold
def manifold_dp(vol_flow):
    vol_flow = vol_flow.to(u.gallon / u.minute)
    
    #10-4-8T
    #dp = vol_flow.magnitude * 3.33 * u.psi

    dp = vol_flow.magnitude * 2.4 * u.psi
    return dp


# Pressure drop functions -> adopted and slightly modified by pressure drop code
def straight_dp(line_vel, length, d_outer, thic, ar):
    if line_vel.magnitude >= 1:
        d_inner = (d_outer - 2 * thic).to(u.meter)  #inner diameter

        rel_rough = core.relative_roughness(d_inner.magnitude, ar.to(u.meter).magnitude)   #relative roughness

        re_eth = core.Reynolds(D=d_inner.magnitude, V=line_vel.to(u.meter / u.second).magnitude, nu=kinetic_mu.magnitude)
        fric_eth = fittings.friction_factor(re_eth, eD=rel_rough)

        # Calcualation of loss coefficient
        K = core.K_from_f(fric_eth, length.to(u.meter).magnitude, d_inner.to(u.meter).magnitude)
        dp = core.dP_from_K(K, rho=density.to(u.kilogram / u.meter**3).magnitude, V=line_vel.to(u.meter / u.second).magnitude) * u.pascal
        dp = dp.to(u.psi)
    else: 
        dp = 0 * u.psi
    return dp

def bend_dp(line_vel, angle, rad, d_outer, thic, ar):
    if line_vel.magnitude >= 1:

        d_inner = (d_outer - 2 * thic).to(u.meter)  #inner diameter

        rel_rough = core.relative_roughness(d_inner.magnitude, ar.to(u.meter).magnitude)   #relative roughness

        re = core.Reynolds(D=d_inner.magnitude, V=line_vel.to(u.meter / u.second).magnitude, nu=kinetic_mu.magnitude)
        fric = fittings.friction_factor(re, eD=rel_rough)

        # Calcualation of loss coefficient
        K = fittings.bend_rounded(
            Di=d_inner.magnitude,  # Inner diameter
            angle=angle.to(u.degree).magnitude,  # degree
            fd=fric, 
            rc=rad.to(u.meter).magnitude,  # Use magnitude of radius in meters
            Re=re, 
            method='Crane'
        )
        dp = core.dP_from_K(K, rho=density.to(u.kilogram / (u.meter**3)).magnitude, V=line_vel.to(u.meter / u.second).magnitude) * u.pascal
        dp = dp.to(u.psi)
    else:
        dp = 0 * u.psi
    return dp

def valve_dp(line_vel, cv, d_valve):
    #if line_vel.magnitude >= 1:
        d_valve = d_valve.to('meter')

        # Calcualation of loss coefficient
        K = fittings.Cv_to_K(Cv = cv, D = d_valve.magnitude)
        dp = core.dP_from_K(K, rho=density.to(u.kilogram / (u.meter**3)).magnitude, V=line_vel.to(u.meter / u.second).magnitude) * u.pascal
        dp = dp.to(u.psi)
    # else:
    #     dp = 0 * u.psi
        return dp


if __name__ == "__main__":
    main()