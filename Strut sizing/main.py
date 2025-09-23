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
    drag = 0.5*(1.225*(velocity**2)/2)*area
    return drag


D = dragCalc(maxV) # proportion of the drag force above a certain point
m = 2

T = -m*(a+g)-D
print(f"Estimated compressive load per strut (uppers) = {T/3}N")
