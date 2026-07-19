import math
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
Cv = 0.8
P1 = 3300
T = 300 * c.KELVIN2RANK
Gs = 0.967

qdot = Cv * 13.61*P1*math.sqrt(1/(T*Gs))
print(qdot)