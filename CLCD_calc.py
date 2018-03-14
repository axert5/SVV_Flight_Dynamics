#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 14:30:31 2018

@author: MatsKlijn
"""

from numpy import *
import matplotlib.pyplot as plt
from WeightBalance import cg

S = 30.00                                                                      #m^2
lbs_to_kg = 0.45359237
g = 9.81                                                                       #gravity constant
fuel_start = 4100*lbs_to_kg*g                                                  #fuel at start in N
weight_plane_empty = cg(0,0)[0]*g                                                    #Weight empty aircraft in N

weights_passengers_notfloat = [82, 92, 60, 60, 77, 69, 67, 95, 84]             #in kg

weights_passengers = []
for i in range(len(weights_passengers_notfloat)):
    weights_passengers_notfloat[i] = float(weights_passengers_notfloat[i])
    weights_passengers.append(weights_passengers_notfloat[i]*g)                #in N
    
#print('weights passengers', weights_passengers)

fuel_used_lbs = [385, 416, 442, 471, 502, 544]                                 #fuel used in lbs

fuel_used = []                                                                 #fuel used in kg

for i in range(len(fuel_used_lbs)):
    fuel_used_lbs[i] = float(fuel_used_lbs[i])
    fuel_used.append((fuel_used_lbs[i]*lbs_to_kg)*g)                           #in N

#print('fuel used', fuel_used)

fuel_stationarystart = []                                                      #Total fuel in aircraft at start of stationary test i

for i in range(len(fuel_used)):
    fuel_diff = fuel_start - fuel_used[i]
    fuel_stationarystart.append(fuel_diff)
    
#print('fuel stationary', fuel_stationary)

lift_stationary = []                                                           #Lift at start of stationary test i
for i in range(len(fuel_stationarystart)):
    lift_stationary.append(fuel_stationarystart[i] + weight_plane_empty)
    
#print('lift value', lift_stationary)


V_kts = [248, 221, 191, 161, 141, 117]                                         #kts
kts_to_ms = 0.51444444444                                                      #kts to m/s
V = []

for i in range(len(V_kts)):
    V.append(V_kts[i]*kts_to_ms)

#print (V)
    
temp_0 = 273.15                                                                #standard temperature Kelvin at sea level
total_temp_stationary = [7.3, 2.5, 1.0, -0.2, -1.0, -2.5]                      #Celsius of temperatures during test measured

for i in range(len(total_temp_stationary)):
    total_temp_stationary[i] = array(total_temp_stationary[i] + temp_0)        #Temperature in Kelvin per  test
    
#print(total_temp_stationary)
    
### Rho values from Nick's program

from rho import rho1
hp = [2133.6, 2133.6, 2133.6, 2133.6, 2130.552, 2097.024]
Vc = [126.0966667, 111.6922222, 96.25888888, 80.82555555, 70.53666666, 58.18999999]
Ttot = [278.15, 275.65, 274.15, 272.95, 272.15, 270.65]
rho = []

for i in range(len(hp)):
    rho2 = rho1(hp[i],Vc[i],Ttot[i])
    rho.append(rho2)
print(rho)

### C_L Calculation
    
C_L = []

for i in range(len(lift_stationary)):
    C_L_formula = lift_stationary[i]/(0.5*rho[i]*V[i]**2*S)                    #C_L formula
    C_L.append(C_L_formula)                                                    #C_L values for stationary flight data


#thrust_left = [3695.85, 2842.33, 2483.98, 2097.37, 1651.94, 2037.49]
#thrust_right = [3688.72, 2919.02, 2612.52, 2217.44, 1818.36, 2236.48]
thrust_total = [7236.03, 5761.35, 5096.5, 4314.81, 3470.3, 4273.97]

#for i in range(len(thrust_left)):
#    thrust_total.append(thrust_left[i]+thrust_right[i])

#print (thrust_total)

C_D = []

for i in range(len(thrust_total)):
    C_D_formula = thrust_total[i]/(0.5*rho[i]*V[i]**2*S)                       #C_D formula
    C_D.append(C_D_formula)                                                    #C_D values for stationary flight data

print(C_D)

### Alpha values 

alpha_rad = array([0.02792526803, 0.04188790205, 0.06108652382, 0.09424777961, 0.1291543646, 0.193731547])

plt.plot(C_D, C_L, 'ro')
plt.xlabel('C_D')
plt.ylabel('C_L')
plt.show()

plt.plot(alpha_rad, C_L, 'b')
plt.xlabel('alpha')
plt.ylabel('C_L')
plt.show()

plt.plot(alpha_rad, C_D, 'g')
plt.xlabel('alpha')
plt.ylabel('C_D')
plt.show()

    
    
  







  