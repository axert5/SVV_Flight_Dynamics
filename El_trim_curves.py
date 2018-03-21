#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 09:34:38 2018

@author: MatsKlijn
"""

g = 9.81
rho0 = 1.225    
Ws = 60500                                                      #Standard aircraft weight (N)
W_empty = 3364                                                  #kg
W_fuel_start = 4100*0.45359237                                  #kg
gamma = 1.4
R = 287.05
diameter_jet = 0.686                                             #m
S = 30.00                                                       #m^2
p0 = 101325     
lamda = -0.0065
T0 = 288.15                                            

from numpy import *
from CLCD_calc import V_tas, lift_stationary
from Cit_par import Cm0, Cma, CmTc, Cmde
from rho import rho1
from matplotlib import pyplot as plt

#DATA
W_passengers = sum(array([82, 92, 60, 60, 77, 69, 67, 95, 84]))             #Weight in N
hp = array([7030, 7270, 7500, 7740, 7020, 6680, 6050])*0.3048                  #altitude in m
V_ias = array([165, 155, 145, 136, 175, 186, 196])*0.51444444444               #m/s
alpha = array([5.1, 5.9, 6.9, 7.9, 4.4, 3.8, 3.3])                             #deg
delta_e = array([0, -0.4, -0.8, -1.3, 0.4, 0.7, 1])                            #deg
delta_e_trim = array([2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8])                      #deg
stick_force = array([0, -21, -26, -44, 26, 55, 78])                            #N
fuel_flow_left = array([466, 462, 458, 453, 470, 479, 489])*0.000125997881     #kg/s
fuel_flow_right = array([490, 486, 483,478, 493, 500, 510])*0.000125997881     #kg/s
fuel_used = array([632, 646, 666, 686, 707, 726, 760])*0.45359237              #kg
total_air_temp = array([-0.2, -1.2, -2.2, -3, 0.2, 1.5, 2.8]) + 273.15         #Kelvin
thrust_left = array([2007.77, 2045.94, 2086.33, 2112.19, 1982.68, 1971.24, 1960.71])
thrust_right = array([2185.43, 2226.46, 2275.74, 2306.23, 2149.08, 2121.65, 2108.01])
thrust_total = array(thrust_left+thrust_right)
thrust_standard = array([2786.86, 2899.84, 3018.24, 3126.6, 2699.72, 2575.56, 2442])
delta_e = array([0, -0.4, -0.8, -1.3, 0.4, 0.7, 1])*(math.pi/180.)

#Weights at start tests
W_start = []
for i in range(len(fuel_used)):
    W_start.append((W_empty + W_passengers + W_fuel_start - fuel_used[i])*g)


#Equivalent Airspeed (V_e)
Vc = []
for i in range(len(V_ias)):
    Vc.append(V_ias[i] - 2*0.51444444444)
    
###Pressure calculation
p = []

for i in range(len(hp)):
    p1 = p0*(1 + (lamda*hp[i])/T0)**(-g/(R*lamda))
    p.append(p1)

#Mach  
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

#Calibrated Temp
Temp_c = []
for i in range(len(M)):
    Temp_c.append(total_air_temp[i]/(1+((gamma-1)/2)*M[i]**2))

#Angle of Attack
aoa = []
for i in range(len(M)):
    aoa.append(sqrt(gamma*R*Temp_c[i]))

#True Airspeed
V_tas = []
for i in range(len(M)):
     V_tas.append(aoa[i]*M[i])

rho = []
for i in range(len(hp)):
    rho.append(rho1(hp[i],Vc[i],total_air_temp[i]))

#Reduced Equivalent Airspeed
V_e = []
for i in range(len(rho)):
    V_e.append(V_tas[i]*(sqrt(rho[i])/rho0)*(Ws/W_start[i]))
    
#Dimensionless Thrust Coefficient (Tc)
thrust_c = []
for i in range(len(thrust_total)):
    thrust_c.append(thrust_total[i]/(0.5*rho[i]*V_tas[i]**2*diameter_jet**2)) 

#Standard thrust measurements

thrust_cs = []
for i in range(len(thrust_standard)):
    thrust_cs.append(thrust_standard[i]/(0.5*rho[i]*V_tas[i]**2*diameter_jet**2))
    
#Delta_eq Calculation 
for i in range(len(thrust_c)):
    delta_eq = delta_e-(CmTc/Cmde)*(thrust_cs[i]-thrust_c[i])

delta_eqsorted = sort(delta_eq)
V_esorted = sort(V_e)

plt.plot(V_esorted,delta_eqsorted, 'xr')
#plt.plot(V_esorted,delta_eqsorted, 'b')
plt.xlabel('V_e')
plt.ylabel('delta_eq')
z = polyfit(V_esorted, delta_eqsorted, 1)
p = poly1d(z)
plt.plot(V_esorted,p(V_esorted),"black")
plt.grid()
plt.show()

stick_force_star = []
for i in range(len(stick_force)):
    stick_force_star.append(stick_force[i]*(Ws/W_start[i]))

stick_force_star_sorted = sort(stick_force_star)

plt.plot(V_esorted,stick_force_star_sorted, 'xr')
#plt.plot(V_esorted,stick_forcesorted, 'b')
plt.xlabel('V_e')
plt.ylabel('stick_force')
z = polyfit(V_esorted, stick_force_star_sorted, 1)
p = poly1d(z)
plt.plot(V_esorted,p(V_esorted),"black")
plt.grid()
plt.show()

print('Trendline coefficients [a,b] for delta_eq = a*V_e+b', polyfit(V_esorted,delta_eqsorted, 1))
print('Trendline coefficients [a,b] for stick_force = a*V_e+b', polyfit(V_esorted,stick_force_star_sorted, 1))








