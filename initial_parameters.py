# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 11:09:57 2018

@author: Max
"""

import numpy as np
import matplotlib.pyplot as plt

text_file=open("flightdata.txt","r")

#Time, altitude, TAS, Angle of attack, Pitch, mass used left engine, mass used right engine
file = np.loadtxt('flightdata.txt', delimiter = ',', skiprows = 1, usecols = (47, 36, 41, 0, 21, 13, 14 ))
text_file.close()

file[:,2] *= 0.514444
file[:,[3,4]]*= np.pi/180 
file[:,1] *= 0.3048 
file[:,[5,6]] *= 0.45359237 


# Begin time for each maneouvre IN SECONDS
Tshort = 48*60+11
Tphugoid = 52*60 + 42
Tdutch = 49*60 + 33
Taperiodic= 51*60 + 30
Tspiral= 57*60 + 15

times = [Tshort, Tphugoid, Tdutch, Taperiodic, Tspiral]
times_names = ["Short", "Phugoid", "Dutch", "Aperiodic", "Spiral"]


for t in times:
    
    index = np.where(file[:,0] == t)[0]
    
    print ("\n \n \n ----------", t, "---------- \n")
    print (" Start of maneoucre",        file[index, 0], "sec" )
    print ("\n Altitude",                file[index, 1],  "m"  )
    print ("\n TAS",                     file[index, 2],  "m/s")
    print ("\n Angle of attack",         file[index, 3],  "deg")
    print ("\n Pitch angle",             file[index, 4],  "deg")
    print ("\n Mass used left engine",   file[index, 5],  "kg" )
    print ("\n Mass used right engine",  file[index, 6],  "kg" )


