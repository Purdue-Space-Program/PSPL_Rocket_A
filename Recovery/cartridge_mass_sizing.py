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

def main():
    f_nc = calcNoseconeForce(fPins, safetyFactor)
    p_nc = calcNoseconePressure(f_nc, aNosecone)

    vExp_in3 = calcExpandedVolume(vBay, p_nc, p_atm, gamma)
    vExp_L = vExp_in3 / 61.024

    molesCO2 = calcRequiredMoles(p_atmAbs, vExp_L, rGasConstant, tempK)
    massCO2 = calcRequiredMass(molesCO2, molarMassCO2)

    print(f"\nBay Volume: {vBay:.2f} in^3")
    print(f"Nosecone Area: {aNosecone:.2f} in^2")
    print(f"Required Pin Shear Force: {fPins:.2f} lbf")

    print(f"\nNosecone Force (SF={safetyFactor}): {f_nc:.2f} lbf")
    print(f"Nosecone Pressure: {p_nc:.2f} psi")

    print(f"\nExpanded CO2 Volume: {vExp_in3:.2f} in^3")
    print(f"Expanded CO2 Volume: {vExp_L:.4f} L")

    print(f"\nRequired CO2 Moles: {molesCO2:.4f} mol")
    print(f"Required CO2 Mass:  {massCO2:.2f} g")

main()