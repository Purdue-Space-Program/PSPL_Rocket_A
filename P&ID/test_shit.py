import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import vehicle_parameters as vehicle

ISP = vehicle.parameters.ISP
print(f"ISP: {ISP}")
# --> 162.21

mass_distribution = vehicle.wet_mass_distribution
for mass in mass_distribution:
    print(mass)