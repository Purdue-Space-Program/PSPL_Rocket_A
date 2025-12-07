def sheerStrength(material, diameter):
    return (material)*(3.14159*diameter*diameter*0.25)

def bearingStrength(material, diameter, thickness):
    return (material)*(diameter)*(thickness)

def main():
    fSU_SS_316_U = 66000.0
    boltHoleDiameter = 0.25
    boltSheerStrength = sheerStrength(fSU_SS_316_U, boltHoleDiameter)
    print(boltSheerStrength)

    fBRU_A_6061_U = 88000.0
    fBRU_A_6061_Y = 58000.0
    plateThickness = 0.125
    boltBearingStrengthUltimate = bearingStrength(fBRU_A_6061_U, boltHoleDiameter, plateThickness)
    boltBearingStrengthYield = bearingStrength(fBRU_A_6061_Y, boltHoleDiameter, plateThickness)
    #print(f"boltBearingStrengthYield: {boltBearingStrengthYield}")
    print(boltBearingStrengthUltimate)
    print(boltBearingStrengthYield)

main()