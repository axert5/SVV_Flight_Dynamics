# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 11:09:57 2018

@author: Max
"""

import numpy as np
import matplotlib.pyplot as plt
#from Symmetric_flight import symmetric
from control.matlab import *

#statespace = symmetric()
#print (statespace)
#Cols:
#Time, Angle of Attack, Mach, TAS, Pitch, Pitch rate, Angle of roll, Roll rate, Yaw rate, Aileron deflection, Elevator deflection, Rudder deflection, altitude
a = np.loadtxt('flightdata.txt', delimiter = ',', skiprows = 1, usecols = (47,0,39,41,21,26,20,25,27,15,16,17,36,13,14))
a[:,3] *= 0.514444
a[:,[1,4,5,6,7,8,9,10,11]]*= np.pi/180
a[:,12] *= 0.3048
a[:,[13,14]] *= 0.45359237

def Short_period():
    #Short period
    Tshort = 48*60 + 11
    index = np.where(a[:,0]==Tshort)
    short = a[index[0][0]:(index[0][0]+(13*10)),:]
    short[:,0] = short[:,0]-short[0,0]
    
    
    maximum = np.where(short[:,4]==max(short[:,4]))
    minimum = np.where(short[:,4]==min(short[:,4]))
    #print (maximum, minimum)
    Period = 2 * (abs(short[minimum[0][0],0]-short[maximum[0][0],0]))
    
    return short[:,0], short[:,1], short[:,3], short[:,4], short[:,5], short[:,10], short[:,12], short[:,13] , short[:,14]
#print (Period)

#plt.plot(b[:,0],b[:,3])
#plt.plot(b[:,0],b[:,1])
#plt.show()
#plt.plot(b[:,0],b[:,-2])
#plt.show()



#Phugoid
def Phugoid():
    Tphugoid = 52*60 + 42
    index = np.where(a[:,0]==Tphugoid)
    b = a[index[0][0]:(index[0][0]+(155*10)),:]
    b[:,0] = b[:,0]-b[0,0]
    
    maximumtime = np.where(b[250:,4]==max(b[250:,4]))
    minimumtime = np.where(b[250:,4]==min(b[250:,4]))

    average = np.average(b[:,4])
    
    
    maximum = max(b[250:,4])
    minimum = min(b[250:,4])
    differencefactor = abs(maximum-average)/abs(minimum-average)
    differencetime = abs(minimumtime[0][0]-maximumtime[0][0])/10
    coefficient = np.log(differencefactor)/(differencetime)
    halftime = np.log(0.5)/coefficient
    #print (halftime)
    
    Period = 2 * (abs(b[minimumtime[0][0],0]-(b[maximumtime[0][0],0])))
    return b[:,0], b[:,1], b[:,3], b[:,4], b[:,5], b[:,10], b[:,12], b[:,13] , b[:,14]
#print (Period)

#y = lsim(statespace,b[:,-2],b[:,0])

#plt.plot(b[:,0],b[:,4])
#plt.plot(b[:,0],b[:,5])
#plt.plot(b[:,0], y[0][:,3])
#plt.plot(b[:,0], y[0][:,2])
plt.show()
#plt.plot(b[:,0],b[:,-2])
plt.show()


plt.show()

#plt.show()



#Dutch roll
Tdutch = 49*60 + 33
index = np.where(a[:,0]==Tdutch)
b = a[index[0][0]:(index[0][0]+(20*10)),:]
b[:,0] = b[:,0]-b[0,0]

fit = (np.polyfit(b[:,0],b[:,6],1))

#plt.plot(b[:,0],b[:,6])
#plt.plot(b[:,0],b[:,7])
#plt.plot(b[:,0],b[:,8])
#plt.show()
#plt.plot(b[:,0],b[:,-1])
#plt.show()

maximumtime = np.where(b[62:,7]==max(b[62:,7]))
minimumtime = np.where(b[62:,7]==min(b[62:,7]))
#print (np.where(np.logical_and(b[:,7]>=-0.4, b[:,7]<=0.4)))
np.where(np.logical_and(b[:,7]>=-0.5, b[:,7]<=0.5))
#print (maximum, minimum)
average = np.average(b[62:,7])

maximum = max(b[62:,7])
minimum = min(b[62:,7])
differencefactor = abs(maximum-average)/abs(minimum-average)
if differencefactor > 1:
    differencefactor = 1/differencefactor
differencetime = abs(minimumtime[0][0]-maximumtime[0][0])/10
coefficient = np.log(differencefactor)/(differencetime)
halftime = np.log(0.5)/coefficient
#print (halftime)

Period = 2 * (abs(b[minimumtime[0][0],0]-(b[maximumtime[0][0],0])))
#print (Period)

#Aperiodic roll
#Taperiodic= 51*60 + 30
#index = np.where(a[:,0]==Taperiodic)
#b = a[index[0][0]:(index[0][0]+(40*10)),:]
#b[:,0] = b[:,0]-b[0,0]
#
#plt.plot(b[:,0],b[:,6])
#plt.plot(b[:,0],b[:,7])
#plt.plot(b[:,0],b[:,8])
#plt.show()
#plt.plot(b[:,0],b[:,-3])
#plt.show()

#Spiral
#Tspiral= 57*60 + 15
#index = np.where(a[:,0]==Tspiral)
#b = a[index[0][0]:(index[0][0]+(190*10)),:]
#b[:,0] = b[:,0]-b[0,0]
#
#plt.plot(b[:,0],b[:,6])
#plt.plot(b[:,0],b[:,7])
#plt.plot(b[:,0],b[:,8])
#plt.show()
#plt.plot(b[:,0],b[:,-1])
#plt.plot(b[:,0],b[:,-3])
#plt.show()


