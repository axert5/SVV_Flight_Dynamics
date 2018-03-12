#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 14:30:31 2018

@author: MatsKlijn
"""

from numpy import *

S = 30.00                                       #m^2
lbs_to_kg = 0.45359237
g = 9.80665                                     #gravity constant
fuel_start = 4100*lbs_to_kg*g                   #fuel at start in N
weight_plane_empty = 3655*g                     #Weight empty aircraft in N

weights_passengers_notfloat = [82, 92, 60, 60, 77, 69, 67, 95, 84]          #in kg

weights_passengers = []
for i in range(len(weights_passengers_notfloat)):
    weights_passengers_notfloat[i] = float(weights_passengers_notfloat[i])
    weights_passengers.append(weights_passengers_notfloat[i]*g)             #in N
    
#print('weights passengers', weights_passengers)

fuel_used_lbs = [386, 416, 442, 471, 502, 544]                              #fuel used in lbs

fuel_used = []                                                              #fuel used in kg

for i in range(len(fuel_used_lbs)):
    fuel_used_lbs[i] = float(fuel_used_lbs[i])
    fuel_used.append((fuel_used_lbs[i]*lbs_to_kg)*g)                        #in N

#print('fuel used', fuel_used)

fuel_stationarystart = []                       #Total fuel in aircraft at start of stationary test i

for i in range(len(fuel_used)):
    fuel_diff = fuel_start - fuel_used[i]
    fuel_stationarystart.append(fuel_diff)
    
#print('fuel stationary', fuel_stationary)

lift_stationary = []                                        #Lift at start of stationary test i
for i in range(len(fuel_stationarystart)):
    lift_stationary.append(sum(weights_passengers) + fuel_stationarystart[i] + weight_plane_empty)
    
#print('lift value', lift_stationary)


V_kts = [248, 221, 191, 161, 141, 117]                       #kts
kts_to_ms = 0.51444444444                                    #kts to m/s
V = []

for i in range(len(V_kts)):
    V.append(V_kts[i]*kts_to_ms)

#print (V)
    
temp_0 = 273.15                                                 #standard temperature Kelvin at sea level
total_temp_stationary = [7.3, 2.5, 1.0, -0.2, -1.0, -2.5]                #Celsius of temperatures during test measured

for i in range(len(total_temp_stationary)):
    total_temp_stationary[i] = total_temp_stationary[i] + temp_0        #Temperature in Kelvin per  test
    
print(total_temp_stationary)
    
rho_standard = 1.225
temp_standard = 288.15

rho_stationary = []
gamma = 1.4 

for i in range(len(total_temp_stationary)):
    
    rho_stationary.append(rho)

print(rho_stationary)







    
    
    
  







  