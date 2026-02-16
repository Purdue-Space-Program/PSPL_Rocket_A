import numpy as np
import matplotlib.pyplot as plt
import math

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import constants as c
import vehicle_parameters as vehicle
import Chambercontour as chamber


xchamb, ychamb = chamber.nozzle_contour(vehicle.parameters.chamber_throat_diameter, chamber.expansion_ratio, chamber.Lstar, vehicle.parameters.contraction_ratio, chamber.con_angle, vehicle.parameters.chamber_inner_diameter, chamber.filename)

#hoopstress calc

x = np.loadtxt('chamber_contour_inches.csv', delimiter = ',', usecols = 2)
innerr = np.loadtxt('chamber_contour_inches.csv', delimiter = ',', usecols = 0)
sf = 1.5

pf = 2
sof = sf*pf
mechs = 72000 #(psi) 

maxstress = mechs/sf # calculates hoopstress
pressure = 250*pf #psi double nominal chamber pressure
outerr = innerr[0] * math.sqrt((maxstress + pressure)/(maxstress - pressure))
chamberthickness = outerr - innerr[0]

print(f"To ensure a safety factor of 1.5, the thickness of the chamber should be at least {chamberthickness:0.3f} in")

minthickness = []
hoopstressar = []
MoSar = []
realhoopstressar = []
for r in innerr:
    outerr = r + chamberthickness
    hoopstress = ((r**2 * pressure)/(outerr**2 - r**2))*(1+(outerr**2/r**2)) #hoopstress with singular minimum thickness
    hoopstressar.append(hoopstress)

    outerr = r * math.sqrt((maxstress + pressure)/(maxstress - pressure)) #changing thickness for each radius
    minthickness.append((outerr-r))

for r in innerr:
    outerr = r + 0.25
    hoopstress = ((r**2 * pressure)/(outerr**2 - r**2))*(1+(outerr**2/r**2))
    realhoopstressar.append(hoopstress)
    limit_load = hoopstress
    design_load = limit_load*sf*pf   #vehicle.parameters.yield_FOS
    MoS = (mechs/design_load) - 1
    MoSar.append(MoS)

plt.figure(1)
plt.plot(-xchamb, ychamb)
plt.plot(-x, minthickness+ychamb)
plt.xlabel("Axial Distance(in)")
plt.ylabel("Minimum Chamber Thickness Difference(in)")
plt.title("Chamber")
plt.show()

plt.figure(2)
plt.plot(-x, realhoopstressar)
plt.xlabel("Axial Distance(in)")
plt.ylabel("Hoop Stress(psi)")
plt.title("Hoopstress along the chamber" )
plt.show()


plt.figure(3)
plt.plot(-x, minthickness)
plt.xlabel("Axial Distance(in)")
plt.ylabel("Minimum Chamber Thickness(in)")
plt.title("Chamber")
plt.show()

plt.figure(4)
plt.plot(-x, MoSar)
plt.xlabel("Axial Distance(in)")
plt.ylabel("Margin of Safety to Ultimate")
plt.title("Ultimate MoS")
plt.show()
