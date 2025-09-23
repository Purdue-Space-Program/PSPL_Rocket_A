# William Wang
# T(x)=−m(x)(a+g)−D′(x) is the equation I'm trying to use for compressive stresses at each point in the rocket
# this is based on https://wikis.mit.edu/confluence/display/RocketTeam/Calculating+Loads+on+the+Rocket
# units in metric, convert to imperial later

# pathfinder is roughly 9ft tall
mt = 39.0089 # total mass
g = 9.80665 # gravitational accel (m/s^2)
a = 5.9*g # instantaneous acceleration (max, m/s^2)
maxV = 274.4 # .8 mach (max, m/s)


def dragCalc(velocity):
    Cd = 0.55 #assuming 30 degree cone
    area = 0.07075853075 # total surface area of that cone (not including base)
    drag = Cd*(1.225*(velocity**2)/2)*area
    return drag


D = dragCalc(maxV) # proportion of the drag force above a certain point
mtop = 2
mbottom = 30
mmiddle = 15
def loadCalc(velocity, massProportion):
    T = -massProportion * (a + g) - D
    return T
loadBottom = loadCalc(maxV, mbottom)
loadTop = loadCalc(maxV, mtop)
loadMiddle = loadCalc(maxV, mmiddle)


print(f"Estimated compressive load per strut (uppers) = {loadTop/3}N")
print(f"Estimated compressive load per strut (inter-tanks) = {loadMiddle/3}N")
print(f"Estimated compressive load per strut (lowers) = {loadBottom/3}N")
