import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import vehicle_parameters as vehicle

ISP = vehicle.parameters.ISP
print(ISP)
# --> 162.21