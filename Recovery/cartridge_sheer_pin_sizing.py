# import numpy as np

def calculatePressure(mass_g, initialTemperature_K, deltaTemperature_K, radius_in, length_in):
    molesMassCO2_g = 44.01
    moles = mass_g/molesMassCO2_g
    R = 0.0821
    temperature_K = initialTemperature_K + deltaTemperature_K
    volume_L = (3.14159*radius_in*radius_in)*(length_in)*(0.0163871)
    return (((moles)*(R)*(temperature_K))/(volume_L))*14.6959

def calculateUpwardForce(pressure_psi, radius_in, length_in):
    return (pressure_psi)*(3.14159*radius_in*radius_in)

def calculateFrictionForce(pressure_psi, radius_in, length_in, coefficient):
    return (pressure_psi)*(3.14159*radius_in*radius_in)*coefficient

def calculateNetForce(positiveForce, negativeForce):
    return positiveForce - negativeForce

#def sizeSheerPins():
    #return

def main():
    cartridgeMass_g = 26
    roomTemperature_K = 298
    dropTemperature_K = -50
    recoveryBayRadius_in = 5.75/2
    recoveryBayLength_in = 10
    recoveryBayPressure_psi = calculatePressure(cartridgeMass_g, roomTemperature_K, dropTemperature_K, recoveryBayRadius_in, recoveryBayLength_in)
    
    print(f"recoveryBayPressure_psi: {recoveryBayPressure_psi:.2f}")
    
    parachuteRadius_in = 2.45
    upwardPressureForce_lbs = calculateUpwardForce(recoveryBayPressure_psi, parachuteRadius_in, recoveryBayLength_in)
    print(f"upwardPressureForce_lbs: {upwardPressureForce_lbs:.2f}")

    coefficientOfFriction = 0.2
    frictionForce_lbs = calculateFrictionForce(recoveryBayPressure_psi, parachuteRadius_in, recoveryBayLength_in, coefficientOfFriction)
    print(f"frictionForce_lbs: {frictionForce_lbs:.2f}")
    
    netForce_lbs = calculateNetForce(upwardPressureForce_lbs, frictionForce_lbs)
    print(f"netForce_lbs: {netForce_lbs:.2f}")

    #sizeSheerPins()

main()