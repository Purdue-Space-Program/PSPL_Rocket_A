from CoolProp.CoolProp import PropsSI

import math
import numpy as np

def qDotCalc(boilOffRate, density):
    qdot = (boilOffRate/density)*1.2 
    return qdot

def CvCrit(qdot,Gs,T,P1):#from copperhead code
    Cv = qdot/(13.61*P1*math.sqrt(1/(T*Gs)))
    return Cv


#tank parameters
   

P_ox = 250     # psia 
V_ox = 1.95/28.317         # Liters to ft^3
T_ox = 181.782 #deg Rankine

P_fuel = 250     # psia fuel is IPA
V_fuel = 2.55/28.317          # L → m^3
T_fuel = 258.333*1.8 #kelvin to deg Rankine

P_press = 4500  # psia
V_press = 4.7/28.317       # liters to ft^3



#gas properties
NDensity = PropsSI('D','P',101325,'T',294.261,'Nitrogen') #[kg/m^3]
airDensity = PropsSI('D','P',101325,'T',294.261,'Air') #[kg/m^3]
Gs = NDensity/airDensity #specific gravity of nitrogen

loxDensity = PropsSI('D', 'P', 1724000, 'T', 100, 'Oxygen')#[kg/m^3]
loxDensity = loxDensity*0.062428 #[kg/m^3] → lbm/ft^3
boilOffRate_ox = float(input("enter boil off rate for Oxygen in lbm/min: "))
qdot_ox = qDotCalc(boilOffRate_ox, loxDensity)

fuelDensity = 786 #[kg/m^3]
fuelDensity = fuelDensity*0.062428 #[kg/m^3] → lbm/ft^3
boilOffRate_fuel = float(input("enter boil off rate for fuel in lbm/min: "))
qdot_fuel = qDotCalc(boilOffRate_fuel, fuelDensity)

CvCrit_ox = CvCrit(qdot_ox,Gs,T_ox,P_ox)
CvCrit_fuel = CvCrit(qdot_fuel,Gs,T_fuel,P_fuel)

print("Cv Critical for Oxygen: ", CvCrit_ox)
print("Cv Critical for Fuel: ", CvCrit_fuel)

 # loxDensity = PropsSI('D', 'P', 1724000, 'T', 100, 'Oxygen')
  #  waterDensity = PropsSI('D', 'P', 1724000, 'T', 294.261, 'Water')

   # print(loxDensity)
   # print(waterDensity)
   # print(loxDensity/waterDensity)