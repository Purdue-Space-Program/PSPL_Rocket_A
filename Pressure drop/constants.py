import numpy as np
# from CoolProp.CoolProp import PropsSI

# Main function

CONVERGE_TOLERANCE = 0.01  # [kg] Allowable difference between masses for Structures and Propulsion to converge
OUTPUT_PRECISION = 3  # [1] Number of digits to round outputs to

# Conversion Factors

## Pressure

PSI2PA = 6894.76  # [Pa/psi] Conversion factor from psi to Pa
PA2PSI = 1 / PSI2PA  # [psi/Pa] Conversion factor from Pa to psi

ATM2PA = 101325  # [Pa/atm] Conversion factor from atm to Pa
PA2ATM = 1 / ATM2PA  # [atm/Pa] Conversion factor from Pa to atm

PA2BAR = 1e-5  # [bar/Pa] Conversion factor from Pa to bar
BAR2PA = 1 / PA2BAR  # [Pa/bar] Conversion factor from bar to Pa

## Temperature

RANK2KELVIN = 5 / 9  # [K/R] Conversion factor from R to K
KELVIN2RANK = 1 / RANK2KELVIN  # [R/K] Conversion factor from K to R
T_AMBIENT = 290  # [K] Ambient temperature (62 F)

## Mass

LB2KG = 0.453592  # [kg/lbm] Conversion factor from lbm to kg
KG2LB = 1 / LB2KG  # [lbm/kg] Conversion factor from kg to lbm

G2KG = 0.001  # [kg/g] Conversion factor from g to kg
KG2G = 1 / G2KG  # [g/kg] Conversion factor from kg to g

## Length

IN2M = 0.0254  # [m/in] Conversion factor from in to m
M2IN = 1 / IN2M  # [in/m] Conversion factor from m to in

M2FT = 3.28084  # [ft/m] Conversion factor from m to ft
FT2M = 1 / M2FT  # [m/ft] Conversion factor from ft to m

## Area

IN22M2 = IN2M**2  # [m^2/in^2] Conversion factor from in^2 to m^2
M22IN2 = 1 / IN22M2  # [in^2/m^2] Conversion factor from m^2 to in^2

## Volume

L2M3 = 0.001  # [m^3/l] Conversion factor from L to m^3
M32L = 1 / L2M3  # [1/m^3] Conversion factor from m^3 to L

M32IN3 = 61023.7  # [in^3/m^3] Conversion factor from m^3 to in^3
IN32M3 = 1 / M32IN3  # [m^3/in^3] Conversion factor from in^3 to m^3

M32FT3 = 35.3147  # [ft^3/m^3] Conversion factor from m^3 to ft^3
FT32M3 = 1 / M32FT3  # [m^3/ft^3] Conversion factor from ft^3 to m^3

## Force

N2LBF = 0.224809  # [lbf/N] Conversion factor from N to lbf
LBF2N = 1 / N2LBF  # [N/lbf] Conversion factor from lbf to N

# Speed

RPM2RADS = (2 * np.pi) / 60

### Material Properties

# 6000-Series Aluminum (https://asm.matweb.com/search/specificmaterial.asp?bassnum=ma6061t6)

DENSITY_AL = 2700  # [kg/m^3] Density
YIELD_STRENGTH_AL = 276 * 10**6  # [Pa] Yield strength
ULTIMATE_STRENGTH_AL = 310 * 10**6  # [Pa] Ultimate tensile strength
YOUNGS_MODULUS = 68.9 * 10**9  # [Pa] Modulus of elasticity
POISSON_RATIO_AL = 0.33  # [1] Poisson's ratio

# Inconel 718 (https://asm.matweb.com/search/specificmaterial.asp?bassnum=ninc34)

DENSITY_INCO = 8190  # [kg/m^3] Density

## Carbon Fiber

DENSITY_CF = 1790  # [kg/m^3] HexTow AS4 Carbon Fiber Density

## 316 Stainless Steel

DENSITY_SS316 = 7980  # [kg/m^3] Density

# Fluids Constants

VENTURI_DP_RATIO = (
    1 / 1.2
)  # [1] venturi outlet pressure / tank pressure, based on NASA venturi documentation & testing
MISC_DP_RATIO = (
    1 / 1.15
)  # [1] outlet pressure / inlet pressure from misc components, based on CMS numbers

FILM_PERCENT = 10  # [%] Percent of fuel mass flow dedicated to film cooling [WAITING ON ANDREW FOR SOURCE]
RESIDUAL_PERCENT = (
    7  # [%] Percent of propellant mass dedicated to residuals, based on CMS values
)
ULLAGE_PERCENT = (
    10  # [%] Percent of tank volume dedicated to ullage, based on CMS values
)

R_PROP = (
    1 - ULLAGE_PERCENT / 100
)  # [1] Ratio of total tank volume to total propellant volume

BURNOUT_PRESSURE_RATIO = (
    2  # [1] COPV burnout pressure / tank pressure to ensure choked flow
)

K_PRESSURIZATION = 0.64  # [1] Ratio of ideal tank volume to actual tank volume, calculated based on average from past rocket designs
HE_GAS_CONSTANT = 2077.1  # [J/kgK] Helium gas constant
N2_GAS_CONSTANT = 296.80  # [J/kgK] Nitrogen gas constant

FILL_PRESSURE = 60  # [psi] Tank pressure during fill, based on CMS fill ops

# Prop Constants

CHAMBER_WALL_THICKNESS = 0.125 * IN2M  # [in] chamber wall thickness
CHAMBER_FLANGE_WIDTH = 1 * IN2M  # [m] chamber flange thickness

INJECTOR_DP_CHAMBER = (
    0.2  # [1] pressure drop / chamber pressure, based on past rockets & RPE
)
REGEN_DP_CHAMBER = (
    0.5  # [1] pressure drop / chamber pressure, conservatively based on RPE
)

RUNLINE_OD = 0.75 * IN2M
RUNLINE_WALL_THICKNESS = 0.065 * IN2M

# Propellant Properties

WATER_PERCENTAGE = 0  # [1] Percentage of water in the ethanol & IPA mixtures

DENSITY_ETHANOL = 785  # [kg/m^3] Ethanol density at STP (https://www.engineeringtoolbox.com/ethanol-ethyl-alcohol-properties-C2H6O-d_2027.html)
DENSITY_JET_A = 807  # [kg/m^3] Jet-A density at STP (https://www.code7700.com/pdfs/bd700_fuel_density.pdf)
DENSITY_IPA = 785  # [kg/m^3] Isopropyl alcohol density at STP (https://macro.lsu.edu/HowTo/solvents/IPA.htm)
DENSITY_METHANOL = 792  # [kg/m^3] Methanol density at STP (https://www.engineeringtoolbox.com/methanol-d_937.html)
DENSITY_WATER = 1000  # [kg/m^3] Water density at STP (https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html)
DENSITY_GASOLINE = 703  # [kg/m^3] Gasoline density at STP (https://www.engineeringtoolbox.com/gasoline-d_1592.html)
DENSITY_E98 = 794
# DENSITY_ETHANOL = PropsSI("D", "T", 290, "P", 101325, "ethanol")

MAX_POWER = 12000  # max pump power [W]



# Indiana Constants 🧑‍🌾
INDIANA_ALTITUDE = 216 # [m] altitude at ACRE https://ag.purdue.edu/indiana-state-climate/purdue-mesonet/site-metadata/
TRAILER_RAIL_HEIGHT = 20 * FT2M # [m] height of launch trailer rail https://arc.aiaa.org/doi/epdf/10.2514/6.2009-4842 (what a fucking find)

# FAR Constants
FAR_ALTITUDE = 615.09  # [m] altitude of FAR launch site
FAR_RAIL_HEIGHT = 18.29  # [m] height of the rail

# Components

# LIPO Battery Constants

LIPO_CELL_MASS = 0.400  # [kg] mass of a LiPo cell
LIPO_CELL_VOLTAGE = 22.2  # [V] nominal voltage of a LiPo cell
LIPO_CELL_DISCHARGE_CURRENT = 80  # [A] maximum discharge current of a LiPo cell

# Misc
GRAVITY = 9.81  # [m/s^2] acceleration due to gravity

# Fins
FIN_THICKNESS = 0.16 * IN2M  # not sure if this is valid, should discuss further [m]
