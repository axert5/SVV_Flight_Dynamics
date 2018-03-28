# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 15:42:57 2018

@author: Max
"""

import numpy as np
gamma = 1.4
p0 = 101325
rho0 = 1.225
Lambda = -0.0065
g0 = 9.81
R = 287.058
T0 = 288.15
a = np.loadtxt('clcd.csv', delimiter = ',', skiprows = 1, usecols = (3,4,6,7,9))
thrusttable = np.ones((6,5))
thrusttable[:,0] = a[:,0]*0.3048
P = p0*(1 + (thrusttable[:,0]*Lambda)/T0)**(-g0/(R*Lambda))


M = np.sqrt((2/(1.4-1))*((1 + (p0/P)*((1 + (gamma-1)/(2*gamma) * (rho0)/(p0) *((a[:,1]-2)*0.51444444444)**2)**(gamma/(gamma-1))-1))**((gamma-1)/gamma) -1))
thrusttable[:,1] = M


print (thrusttable)