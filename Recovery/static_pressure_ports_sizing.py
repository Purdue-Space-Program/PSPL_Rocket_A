# Sizing Recovery Bay Static Pressure Ports

import math

# Obtained from CAD: Recovery Bay Layout
boreOuterRadius = 1.0156

bottomFaceRadius = 2.1495
bottomFaceHeight = 1.6000

middleFaceRadius = bottomFaceRadius
middleFaceChord = 2.559
middleFaceExt = 0.0777
middleFaceHeight = 0.3960

topFaceRadius = 2.4495
topFaceHeight = 0.2500

def calcCircArea(radius):
    return math.pi * radius**2

def calcCircSegExtArea(radius, chordLength, extLength):
    theta = 2 * math.asin(chordLength / (2 * radius))
    circSegArea = radius**2 * (theta - math.sin(theta))

    extArea = extLength * chordLength

    circSegExtArea = circSegArea + extArea
    return circSegExtArea
    
def calcVolume(area, height):
    return area * height

def calcSinglePortDiam(airVolume):
    if airVolume < 100:
        singlePortDiam = airVolume / 400
    else:
        singlePortDiam = 2 * math.sqrt(airVolume / 6397.71)
    return singlePortDiam

def calcMultiPortDiam(singlePortArea, numHoles):
    return 2 * math.sqrt((singlePortArea / numHoles) / math.pi)

def main():
    numHoles = 3 # Number of ports (may change)

    outer_boreArea = calcCircArea(boreOuterRadius)

    bottomFaceTotArea = calcCircArea(bottomFaceRadius)
    topFaceArea = calcCircArea(topFaceRadius)

    middleFaceAreaTrue = calcCircSegExtArea(middleFaceRadius, middleFaceChord, middleFaceExt)
    bottomFaceAreaTrue = bottomFaceTotArea - (outer_boreArea + middleFaceAreaTrue)
    topFaceAreaTrue = topFaceArea - outer_boreArea

    bottomFaceVol = calcVolume(bottomFaceAreaTrue, bottomFaceHeight)
    middleFaceVol = calcVolume(middleFaceAreaTrue, middleFaceHeight)
    topFaceVol = calcVolume(topFaceAreaTrue, topFaceHeight)

    totalAirVolume = bottomFaceVol + middleFaceVol + topFaceVol
    
    singlePortDiam = calcSinglePortDiam(totalAirVolume)
    singlePortArea = calcCircArea(singlePortDiam / 2)
    multiPortDiam = calcMultiPortDiam(singlePortArea, numHoles)

    print(f"\nBottom Face Area: {bottomFaceAreaTrue:.2f} in^2")
    print(f"Middle Face Area: {middleFaceAreaTrue:.2f} in^2")
    print(f"Top Face Area: {topFaceAreaTrue:.2f} in^2")

    print(f"\nBottom Face Volume: {bottomFaceVol:.2f} in^3")
    print(f"Middle Face Volume: {middleFaceVol:.2f} in^3")
    print(f"Top Face Volume: {topFaceVol:.2f} in^3")

    print(f"\nTotal Air Volume: {totalAirVolume:.2f} in^3")

    print(f"\nSingle Port Diameter: {singlePortDiam:.4f} in")
    print(f"Multi Port Diameter: {multiPortDiam:.4f} in")

main()