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
sof = 1.4
mechs = 72518.9 #(psi) mechanical strength at 500K ~226C

maxstress = mechs/sof # calculates hoopstress
pressure = 500 #psi double nominal chamber pressure
outerr = innerr[0]*math.sqrt((maxstress+pressure)/(maxstress - ((innerr[0]**2 * pressure)/innerr[0]**2)))
chamberthickness = outerr - innerr[0]

print(f"To ensure a safety factor of 1.4, the thickness of the chamber should be at least {chamberthickness*c.M2IN:0.3f} in")

minthickness = []
hoopstressar = []
for i in innerr:
    r = i
    outerr = r + chamberthickness
    hoopstress = ((i**2 * pressure)/(outerr**2 - i**2))*(1+(outerr**2/r**2))
    hoopstressar.append(hoopstress)
    outerr = i*math.sqrt((maxstress+pressure)/(maxstress - ((i**2 * pressure)/i**2)))
    minthickness.append((outerr-i)*c.M2IN)

plt.figure(1)
plt.plot(xchamb, ychamb)
plt.plot(x, minthickness+ychamb)
plt.xlabel("Axial Distance(in)")
plt.ylabel("Minimum Chamber Thickness(in)")
plt.title("Chamber")
plt.show()

plt.figure(2)
plt.plot(x, hoopstressar)
plt.xlabel("Axial Distance(in)")
plt.ylabel("Hoop Stress")
plt.title("Hoopstress along the chamber" )
plt.show()
