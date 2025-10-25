# This code reflects the vehicle parameters as stated on the vehicle parameters page: https://purdue-space-program.atlassian.net/wiki/x/Aw5qYw
# Import this file into your script as a way to easily and reliably pull parameters. 
# Author: David Gustafsson

import constants as c
 
vehicle_parameters = {

    "fuel_name": "isopropyl alcohol", 
    "oxidizer_name": "liquid oxygen",
    
    "chamber_pressure": 150 * c.PSI2PA,      # The target combustion pressure in the engine [Pascals]
    "jet_thrust": 577.52 * c.LBF2N,          # The targeted engine thrust (not accounting for exhaust gas expansion thrust) [Newtons]
    "ISP": 162.21,                           # The estimated ISP of the engine [seconds]
    "mass_flow_rate": 3.56 * c.LB2KG,        # The targeted mass flow rate through the engine [kilograms/second]
    "OF_ratio": 1.0,                         # The target ratio of oxygen to fuel combustion in the engine [dimensionless] 
    "burn_time": 2.24,                       # The estimated burn time of the engine [seconds]
    "contraction_ratio": 3.0,                # The target ratio of chamber area to throat area [dimensionless]
    "exit_pressure": 15.0 * c.PSI2PA,        # The target exit pressure of the exhaust gas [Pascals]
    # "combustion_temperature": 2170,          # The estimated combustion temperature [Kelvin]
    
    "tank_pressure": 250.0 * c.PSI2PA,       # The estimated required tank pressure to sustain the combustion pressure in the engine [Pascals]
    
    "fuel_tank_length": 0.5 * c.FT2M,        # The length of the fuel tank that needs to be filled with fuel (the actual tank may be longer) [meters]
    "fuel_tank_volume": 2.55 * c.L2M3,       # The required volume of fuel needed for the burn time [meter^3]
    "fuel_total_mass": 4.42 * c.LB2KG,       # The required mass of fuel needed for the burn time [kilogram]
    
    "oxidizer_tank_length": 4.57 * c.IN2M,   # The length of the oxidizer tank that needs to be filled with oxidizer (the actual tank may be longer) [meters]
    "oxidizer_tank_volume": 1.95 * c.L2M3,   # The required volume of oxidizer needed for the burn time [meter^3]
    "oxidizer_total_mass": 4.42 * c.LB2KG,   # The required mass of oxidizer needed for the burn time [kilogram]
    
    
    "total_length": 8.80 * c.FT2M,           # The estimated length of the rocket [meter]
    "wet_mass": 77.29 * c.LB2KG,             # The estimated dry mass of the rocket [kilogram]
    "dry_mass": 69.32 * c.LB2KG,             # The estimated dry mass of the rocket [kilogram]
    "estimated_apogee": 3300 * c.FT2M,       # The estimated 1-DOF altitude [meters]
    "off_the_rail_TWR": 7.81,                # The target thrust-to-weight ratio of the rocket off the launch rail [dimensionless]
    "off_the_rail_acceleration": 6.81,       # The target acceleration of the rocket off the launch rail [standard gravity]
    "off_the_rail_velocity": 29.72,          # The target velocity of the rocket off the launch rail [meters/second]
    "max_acceleration": 7.47, # (upwards)    # The maximum acceleration of the rocket during flight [standard gravity]
    # "max_velocity": 0.454,                   # The maximum speed of the rocket during flight [Mach (speed of sound of air)]
    "max_velocity": 155.722,                 # The maximum speed of the rocket during flight [meters/second]
    "total_impulse": 5934,                   # The total impulse of the rocket over the duration of flight [newton seconds]
}