#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 14:30:31 2018

@author: MatsKlijn
"""

from numpy import *
import matplotlib.pyplot as plt
from WeightBalance import cg

S = 30.00  
b = 15.911  
A = b**2/S                                                                  #m^2
lbs_to_kg = 0.45359237
g = 9.81                                                                       #gravity constant
fuel_start = 4100*lbs_to_kg*g                                                  #fuel at start in N
weight_plane_empty = cg(0,0)[0]*g                                                    #Weight empty aircraft in N
R = 287.05
gamma = 1.4
p0 = 101325
T0 = 288.15
lamda = -0.0065
rho0 = 1.225


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


V_kts = [249, 221, 191, 161, 141, 117]                                         #kts
kts_to_ms = 0.51444444444                                                      #kts to m/s
V = []

for i in range(len(V_kts)):
    V.append((V_kts[i]*kts_to_ms))
V=array(V)
Vc = V-2*.51444444444

#print ('Vc is',Vc)
    
temp_0 = 273.15                                                                #standard temperature Kelvin at sea level
total_temp_stationary = [7.3, 2.5, 1.0, -0.2, -1.0, -2.5]                      #Celsius of temperatures during test measured

for i in range(len(total_temp_stationary)):
    total_temp_stationary[i] = array(total_temp_stationary[i] + temp_0)        #Temperature in Kelvin per  test
    
#print(total_temp_stationary)
    
### Rho values from Nick's program

from rho import rho1
hp = [2133.6, 2133.6, 2133.6, 2133.6, 2130.552, 2097.024]
#Ttot = [278.15, 275.65, 274.15, 272.95, 272.15, 270.65]
rho = []

for i in range(len(hp)):
    rho2 = rho1(hp[i],Vc[i],total_temp_stationary[i])
    rho.append(rho2)
#print('rho is',rho)

###Pressure calculation

p = []

for i in range(len(hp)):
    p1 = p0*(1 + (lamda*hp[i])/T0)**(-g/(R*lamda))
    p.append(p1)

#print ('p is', p)
    
#Mach number calculation

M = []

for i in range(len(p)):
    M1 = 1+((gamma-1)/(2*gamma))*(rho0/p0)*Vc[i]**2
    M2 = M1**(gamma/(gamma-1))
    M3 = M2-1
    M4 = 1+(p0/p[i])*M3
    M5 = M4**((gamma-1)/gamma)
    M6 = M5-1
    M7 = (2/(gamma-1))*M6
    M8 = sqrt(M7)
    M.append(M8)
            
#print ('M is', M)

#Calibrated Temp

Temp_c = []

for i in range(len(M)):
    Tc1 = total_temp_stationary[i]/(1+((gamma-1)/2)*M[i]**2)
    Temp_c.append(Tc1)

#Angle of attack
aoa = []
for i in range(len(M)):
    aoa1 = sqrt(gamma*R*Temp_c[i])
    aoa.append(aoa1)
    
#True Airspeed
V_tas = []
for i in range(len(M)):
     V_tas.append(aoa[i]*M[i])

### C_L Calculation
    
C_L = []
C_Lsquared = []

for i in range(len(lift_stationary)):
    C_L_formula = lift_stationary[i]/(0.5*rho[i]*V_tas[i]**2*S)                #C_L formula
    C_L.append(C_L_formula)                                                    #C_L values for stationary flight data
    C_L2 = (C_L_formula)**2
    C_Lsquared.append(C_L2)                                                    #

#thrust_left = [3695.85, 2842.33, 2483.98, 2097.37, 1651.94, 2037.49]
#thrust_right = [3688.72, 2919.02, 2612.52, 2217.44, 1818.36, 2236.48]
thrust_total = [7236.03, 5761.35, 5096.5, 4314.81, 3470.3, 4273.97]

#for i in range(len(thrust_left)):
#    thrust_total.append(thrust_left[i]+thrust_right[i])

#print (thrust_total)

C_D = []

for i in range(len(thrust_total)):
    C_D_formula = thrust_total[i]/(0.5*rho[i]*V_tas[i]**2*S)                   #C_D formula
    C_D.append(C_D_formula)                                                    #C_D values for stationary flight data


#corrected coefficients
corrected_C_D = []
C_D_0 = 0.022
e = 1/(pi*A*0.043)

for i in range(len(C_Lsquared)):
    corrected_C_D_formula = C_D_0 + C_Lsquared[i]/(pi*e*A)                 #C_D formula
    corrected_C_D.append(corrected_C_D_formula)                                                    #C_D values for stationary flight data

### Alpha values 

alpha_rad = array([0.02792526803, 0.04188790205, 0.06108652382, 0.09424777961, 0.1291543646, 0.193731547])

plt.figure(1)

plt.subplot(511)
plt.title('C_L vs C_D ')
plt.plot(corrected_C_D, C_L, 'r')
plt.xlabel('C_D')
plt.ylabel('C_L')

plt.subplot(512)
plt.title('C_L vs alpha')
plt.plot(alpha_rad, C_L, 'b')
plt.xlabel('alpha')
plt.ylabel('C_L')

plt.subplot(513)
plt.title('C_D vs alpha')
plt.plot(alpha_rad, corrected_C_D, 'g')
plt.xlabel('alpha')
plt.ylabel('C_D')

plt.subplot(514)
plt.title('C_L^2 vs C_D')
plt.plot(C_Lsquared,corrected_C_D, 'y')
plt.xlabel('C_L^2')
plt.ylabel('C_D')
#<<<<<<< HEAD
plt.show()
'''
  
=======

plt.subplot(515)
plt.title('C_D vs alpha')
plt.plot(alpha_rad, corrected_C_D, 'g')
plt.xlabel('alpha')
plt.ylabel('C_D')

plt.show()    
    
 


>>>>>>> 3402d9fb1c7a939ddfd89a411fc55cc2616abfd5
'''






  