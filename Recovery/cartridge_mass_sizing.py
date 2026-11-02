# Sizing Recovery Bay CO2 Cartridge Mass

import math

# Obtained from Recovery Bay Layout CAD
rInner = 2.75 # [in]
lBay = 30.0   # [in]

# Recovery Bay Parameters
vBay = math.pi * (rInner**2) * lBay # [cuin]
aNosecone = math.pi * (rInner**2)   # [sqin]

# Shear Pins Parameters
numPins = 3
forcePerPin = 40.0             # [lbf]
fPins = numPins * forcePerPin  # [lbf]

# Deployment & Thermodynamic Constraints
safetyFactor = 1.5    # [N/A]
p_atm = 14.7          # [psi]
p_atmAbs = 1.0        # [atm]
gamma = 1.3           # [N/A]
tempK = 273 - 50      # [K]
rGasConstant = 0.0821 # [L atm / (mol K)]
molarMassCO2 = 44.01  # [g/mol CO2]

# Chosen Cartridge Mass
chosenMass = 24       # [g]

def calcNoseconeForce(shearForce, sf):
    return shearForce * sf

def calcNoseconePressure(force, area):
    return force / area

def calcExpandedVolume(vBay, pNC, pAtm, gamma):
    return vBay * ((pNC / pAtm) ** (1 / gamma))

def calcRequiredMoles(pAbs, vExpL, rConst, temp):
    return (pAbs * vExpL) / (rConst * temp)

def calcRequiredMass(moles, molarMass):
    return moles * molarMass

def calcProvidedPressure(p_atm, vExp, vBay, gamma):
    return p_atm * ((vExp / vBay) ** gamma)

def main():
    f_nc = calcNoseconeForce(fPins, safetyFactor)
    p_nc = calcNoseconePressure(f_nc, aNosecone)

    min_vExp_in3 = calcExpandedVolume(vBay, p_nc, p_atm, gamma)
    min_vExp_L = min_vExp_in3 / 61.024

    molesCO2 = calcRequiredMoles(p_atmAbs, min_vExp_L, rGasConstant, tempK)
    massCO2 = calcRequiredMass(molesCO2, molarMassCO2)

    chosen_vExp_in3 = (chosenMass / massCO2) * min_vExp_in3
    p_provided = calcProvidedPressure(p_atm, chosen_vExp_in3, vBay, gamma)

    print(f"\nBay Volume: {vBay:.2f} in^3")
    print(f"Nosecone Area: {aNosecone:.2f} in^2")
    print(f"Required Pin Shear Force: {fPins:.2f} lbf")

    print(f"\nNosecone Force (SF={safetyFactor}): {f_nc:.2f} lbf")
    print(f"Nosecone Pressure: {p_nc:.2f} psi")

    print(f"\nExpanded CO2 Volume (16.86g): {min_vExp_in3:.2f} in^3")
    print(f"Expanded CO2 Volume (24g): {chosen_vExp_in3:.2f} in3")
    print(f"Expanded CO2 Volume (16.86g: {min_vExp_L:.4f} L")

    print(f"\nMinimum Required CO2 Moles: {molesCO2:.4f} mol")
    print(f"Minimum Required CO2 Mass:  {massCO2:.2f} g")

    print(f"\nProvided Pressure: {p_provided:.2f} psi")

main()