# William Wang
# T(x)=−m(x)(a+g)−D′(x) is the equation I'm trying to use for compressive stresses at each point in the rocket
# this is based on https://wikis.mit.edu/confluence/display/RocketTeam/Calculating+Loads+on+the+Rocket
# units in metric, convert to imperial later
import math
# pathfinder is roughly 9ft tall
mt = 39.0089 # total mass
g = 9.80665 # gravitational accel (m/s^2)
a = 5.9*g # instantaneous acceleration (max, m/s^2)
maxV = 274.4 # .8 mach (max, m/s)
d = 6/39.37007874

def dragCalc(velocity):
    Cd = 0.55 #assuming 30 degree cone
    area = math.pi * (d/2)**2 # total surface area of that cone (not including base)
    drag = Cd*(1.225*(velocity**2)/2)*area
    return drag

mtop = 2
mbottom = 30
mmiddle = 15
def loadCalc(velocity, massProportion):
    D = dragCalc(velocity) # proportion of the drag force above a certain point
    T = -massProportion * (a + g) - D
    return T
loadBottom = loadCalc(maxV, mbottom)
loadTop = loadCalc(maxV, mtop)
loadMiddle = loadCalc(maxV, mmiddle)

if input("use default for mass (39.0089kg)? (y/n)\n> ") == "n":
    mt = float(input("enter mass (kg)\n>"))
if input("use default for max. accel (5.9G)? (y/n)\n> ") == "n":
    a = float(input("enter max acceleration (g)\n>"))
if input("use default for max. velocity (274.4m/s)? (y/n)\n> ") == "n":
    maxV = float(input("enter max velocity (m/s)\n>"))
if input("use default for rocket diameter (6in)? (y/n)\n> ") == "n":
    d = float(input("enter diameter (in)\n>"))/39.37007874


print(f"Estimated compressive load per strut (uppers) = {loadTop/3}N")
print(f"Estimated compressive load per strut (inter-tanks) = {loadMiddle/3}N")
print(f"Estimated compressive load per strut (lowers) = {loadBottom/3}N")
