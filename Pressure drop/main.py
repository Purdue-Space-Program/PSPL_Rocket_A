from fluids import fittings
from fluids import core
from CoolProp import CoolProp as CP
from pint import UnitRegistry
import math

u = UnitRegistry()

def find_drop_ox(K, type, line_velocity_ox):
    dP = core.dP_from_K(K, rho=rho_ox.magnitude, V=line_velocity_ox.magnitude) * u.Pa
    dP.to('psi')
    print(f'LOX LINE: Pressure Drop at {type}: {dP:.2f}')