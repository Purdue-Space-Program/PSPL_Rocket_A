import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import vehicle_parameters
parameters, wet_mass_distribution, _ = vehicle_parameters.main()

ISP = parameters.ISP
print(f"ISP: {ISP}")
# --> 162.21

mass_distribution = wet_mass_distribution
for mass in mass_distribution:
    print(mass)