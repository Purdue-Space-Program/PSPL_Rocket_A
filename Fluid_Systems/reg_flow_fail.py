import math
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
import numpy as np
import matplotlib.pyplot as plt

Cv = 0.8
P1 = 300 * c.BAR2PA * c.PA2PSI
T = 300 * c.KELVIN2RANK
Gs = 0.967

qdot = Cv * 13.61*P1*math.sqrt(1/(T*Gs))
print(qdot)

scfm_history = []
pressure_history = np.linspace(3000, 5000, 100)
for P1 in pressure_history:
    qdot = Cv * 13.61*P1*math.sqrt(1/(T*Gs))
    scfm_history.append(qdot)

plt.plot(pressure_history, scfm_history)
plt.xlabel('Pressure (psia)')
plt.ylabel('Flow Rate (SCFM)')
plt.title('Reg Fail Flow Rate vs Inlet Pressure')
plt.axvline(300 * c.BAR2PA * c.PA2PSI, color='r', linestyle='--', label='COPV MAWP')
plt.axvline(3300, color='b', linestyle='--', label='Nominal Inlet Pressure')
plt.axhline(729.33, color='g', linestyle='--', label='Max Relief Flow Rate at set 475 psig')
plt.axhline(842.38, color='c', linestyle='--', label='Max Relief Flow Rate at set 550 psig')
plt.axhline(729.33 * 3, color='m', linestyle='--', label='3X Max Relief Flow Rate at set 475 psig')
plt.axhline(842.38 * 3, color='y', linestyle='--', label='3X Max Relief Flow Rate at set 550 psig')
plt.ylim(0, 3000)
plt.xlim(3000, 5000)
plt.legend()
plt.show()