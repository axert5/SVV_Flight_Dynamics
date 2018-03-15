#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 11:58:57 2018

@author: MatsKlijn
"""

from numpy import *
from Cit_par import *
from control.matlab import*
import matplotlib.pyplot as plt

C1_asymmetric = matrix([[(CYbdot-2*mub)*b/V0,0,0,0],
                         [0,-b/(2*V0),0,0], 
                         [0,0,(-4*mub*KX2*b**2)/(2*V0**2),(4*mub*KXZ*b**2)/(2*V0**2)],
                         [Cnbdot*b/(V0),0,(4*mub*KXZ*b**2)/(2*V0**2),(-4*mub*KZ2*b**2)/(2*V0**2)]
                         ])

C2_asymmetric = matrix([[CYb, CL, CYp*b/(2*V0),(CYr-4*mub)*b/(2*V0)],
                         [0,0,b/(2*V0),0],
                         [Clb,0,Clp*b/(2*V0),Clr*b/(2*V0)],
                         [Cnb,0,Cnp*b/(2*V0),Cnr*b/(2*V0)]
                         ])

C3_asymmetric = matrix([[CYda,CYdr],
                        [0,0],
                        [Clda,Cldr],
                        [Cnda,Cndr]
                        ])


A_asymmetric = linalg.inv(-C1_asymmetric)*C2_asymmetric
B_asymmetric = linalg.inv(-C1_asymmetric)*C3_asymmetric
C_asymmetric = identity(4)
D_asymmetric = zeros((4,2))

sys_asymmetric=ss(A_asymmetric,B_asymmetric,C_asymmetric,D_asymmetric)


print ('Eigenvalues of A_asymmetric:',linalg.eig(A_asymmetric)[0] )
print ('Eigenvectors of A_asymmetric:',linalg.eig(A_asymmetric)[1] )
#---input to state space-------------------------------------------------------


t=arange(0,10,0.01)

y=impulse(sys_asymmetric,t,input=0)

#----------plotting-----------------------------------------------------------
plt.figure(1)

plt.subplot(411)
plt.title('Sideslip Angle beta (degs)')
plt.plot(t,y[0][:,0]*180/pi,color='c',label='u')

plt.subplot(412)
plt.title('Roll Angle phi (degs)')
plt.plot(t,y[0][:,1]*180/pi, color='r', label='alpha')

plt.subplot(413)
plt.title('Roll Rate p (m/s)')
plt.plot(t,y[0][:,2],color='b',label='theta')

plt.subplot(414)
plt.title('Yaw Rate r (m/s)')
plt.plot(t,y[0][:,3],color='g',label='q')
plt.plot()




