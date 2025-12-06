def calculateUpwardPressureForce(mass_g, volume2_L):
    molesMassN2_g = 28.02
    moles = mass_g/molesMassN2_g
    R = 0.0821
    temperature_K = 298
    pressure1_atm_volume1_L = ((moles)*(R)*(temperature_K))
    pressure2_atm = pressure1_atm_volume1_L/volume2_L
    print(pressure2_atm)
    return (pressure2_atm*14.6959)*(3.14159*2.875*2.875)

def calculateFrictionForce(force, coefficient):
    return coefficient*force

def calculateNetForce(positiveForce, negativeForce):
    return positiveForce - negativeForce

#def sizeSheerPins():
    #return

def main():

    cartridgeMass_g = 40
    recoveryVolume_L = 4.2
    upwardPressureForce_lbs = calculateUpwardPressureForce(cartridgeMass_g, recoveryVolume_L)
    
    coefficientOfFriction = 0.6
    frictionForce_lbs = calculateFrictionForce(upwardPressureForce_lbs, coefficientOfFriction)
    
    netForce_lbs = calculateNetForce(upwardPressureForce_lbs, frictionForce_lbs)
    print(netForce_lbs)

    #sizeSheerPins()

main()