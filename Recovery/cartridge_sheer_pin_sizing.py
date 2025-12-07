def calculateUpwardPressureForce(mass_g, deltaTemperature_K, parachuteRadius_in, radius_in, length_in):
    molesMassCO2_g = 44.01
    moles = mass_g/molesMassCO2_g
    R = 0.0821
    temperature_K = 298-deltaTemperature_K
    volume_L = (3.14159*radius_in*radius_in)*(length_in)*(0.0163871)
    pressure_atm = ((moles)*(R)*(temperature_K))/(volume_L)
    return (pressure_atm*14.6959)*(3.14159*parachuteRadius_in*parachuteRadius_in)

def calculateFrictionForce(force, coefficient):
    return coefficient*force

def calculateNetForce(positiveForce, negativeForce):
    return positiveForce - negativeForce

#def sizeSheerPins():
    #return

def main():

    cartridgeMass_g = 26
    temperatureDrop_K = 50
    parachuteRadius_in = 2.45
    recoveryBayRadius_in = 2.875
    recoveryBayLength_in = 10
    upwardPressureForce_lbs = calculateUpwardPressureForce(cartridgeMass_g, temperatureDrop_K, parachuteRadius_in, recoveryBayRadius_in, recoveryBayLength_in)
    
    coefficientOfFriction = 0.6
    frictionForce_lbs = calculateFrictionForce(upwardPressureForce_lbs, coefficientOfFriction)
    
    netForce_lbs = calculateNetForce(upwardPressureForce_lbs, frictionForce_lbs)
    print(netForce_lbs)

    #sizeSheerPins()

main()