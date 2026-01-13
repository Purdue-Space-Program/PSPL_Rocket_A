molesMassCO2_g = 44.01
R = 0.0821
pi = 3.14159
cuin_to_L = 0.0163871
atm_to_psi = 14.6959
m_to_in = 39.3701
import sys
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from constants import *
import vehicle_parameters as v


def calculatePressure(mass_g, initialTemperature_K, deltaTemperature_K, radius_in, length_in):
    moles = mass_g/molesMassCO2_g
    temperature_K = initialTemperature_K + deltaTemperature_K
    volume_L = (pi*radius_in*radius_in)*(length_in)*(cuin_to_L)
    return (((moles)*(R)*(temperature_K))/(volume_L))*atm_to_psi

def calculateUpwardForce(pressure_psi, radius_in):
    return (pressure_psi)*(pi*radius_in*radius_in)

def calculateFrictionForce(pressure_psi, radius_in, coefficient):
    return (pressure_psi)*(pi*radius_in*radius_in)*coefficient

def calculateNetForce(positiveForce, negativeForce):
    return positiveForce - negativeForce

def calculateSheerPinMOS(calculated, accepted):
    return (calculated - accepted)/accepted

def main():
    cartridgeMass_g = 26
    roomTemperature_K = 273
    dropTemperature_K = -50
    recoveryBayRadius_in = (v.parameters.tube_inner_diameter*m_to_in)/2
    recoveryBayLength_in = v.recovery_bay_length*m_to_in + (v.nosecone_length*m_to_in)/2
    recoveryBayPressure_psi = calculatePressure(cartridgeMass_g, roomTemperature_K, dropTemperature_K, recoveryBayRadius_in, recoveryBayLength_in)
    
    print(f"recoveryBayPressure_psi: {recoveryBayPressure_psi:.2f}")
    
    parachuteRadius_in = 2.45
    upwardPressureForce_lbs = calculateUpwardForce(recoveryBayPressure_psi, parachuteRadius_in)
    print(f"upwardPressureForce_lbs: {upwardPressureForce_lbs:.2f}")

    coefficientOfFriction = 0.2 
    #suspect coefficient of friction
    #no bulkhead nosecone friction?
    frictionForce_lbs = calculateFrictionForce(recoveryBayPressure_psi, parachuteRadius_in, coefficientOfFriction)
    print(f"frictionForce_lbs: {frictionForce_lbs:.2f}")
    
    netForce_lbs = calculateNetForce(upwardPressureForce_lbs, frictionForce_lbs)
    print(f"netForce_lbs: {netForce_lbs:.2f}")
    #figure out what forces shear pins see, and make sure they dont come off when you don't want them too 

    acceptedForce_lbs = 119.2
    MOS = calculateSheerPinMOS(netForce_lbs, acceptedForce_lbs)
    print(f"Margin of Safety: {MOS:.2f}")

main()