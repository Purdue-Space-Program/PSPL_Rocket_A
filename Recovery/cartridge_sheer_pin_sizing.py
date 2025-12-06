def calculateUpwardPressureForce(mass_g, volume2_L):
    molesMassN2_g = 28.02
    moles = mass_g/molesMassN2_g
    R = 0.0821
    temperature_K = 298
    pressure1_atm_volume1_L = ((moles)*(R)*(temperature_K))
    pressure2_atm = pressure1_atm_volume1_L/volume2_L
    return pressure2_atm*(3.14159*2.875*2.875)



def main():

    cartridgeMass_g = 40
    recoveryVolume = 4.2
    upwardPressureForce = calculateUpwardPressureForce(cartridgeMass_g, recoveryVolume)
    
    calculateFrictionForce(recoveryVolume)
    
    calculateNetForce()
    
    sizeSheerPins()

main()