# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 11:11:44 2018

@author: xx
"""
from Cit_par import *
from numpy import*
from control.matlab import*
import matplotlib.pyplot as plt
from scipy.integrate import simps
#State space system

#symmetric case

C1_symmetric= matrix([[-2*muc*c/(V0**2) ,0, 0 ,0],[0,(CZadot-2*muc)*c/V0,0,0],[0,0,-c/V0,0],[0,Cma*c/V0,0,-2*muc*KY2**2*c**2/V0**2]])
C2_symmetric=matrix([[CXu/V0,CXa,CZ0,CXq*c/V0],[CZu/V0,CZa,-CX0,(CZq+2*muc)*c/V0],[0,0,0,c/V0],[Cmu/V0,Cma,0,Cmq*c/V0]])
C3_symmetric=matrix([[CXde],[CZde],[0],[Cmde]])

A_symmetric=linalg.inv(-C1_symmetric)*C2_symmetric
B_symmetric=linalg.inv(-C1_symmetric)*C3_symmetric
C_symmetric=identity(4)
D_symmetric=zeros((4,1))

sys_symmetric=ss(A_symmetric,B_symmetric,C_symmetric,D_symmetric)

#---input to state space-------------------------------------------------------


t=arange(0,10,0.01)
u=[-0.01]*len(t)
y=lsim(sys_symmetric,u,t)

#htemp=y[0][:,0]*sin(y[0][:,2])
#h=simps(htemp,t)

#----------plotting-----------------------------------------------------------
plt.figure(1)

plt.subplot(411)
plt.title('Velocity u (m/s)')
plt.plot(t,y[0][:,0],color='c',label='u')

plt.subplot(412)
plt.title('Angle of Attack alpha (degs)')
plt.plot(t,y[0][:,1]*180/pi, color='r', label='alpha')

plt.subplot(413)
plt.title('Fligth Path Angle theta (degs)')
plt.plot(t,y[0][:,2]*180/pi,color='b',label='theta')

plt.subplot(414)
plt.title('Pitch Rate q (m/s)')
plt.plot(t,y[0][:,3],color='g',label='q')

#plt.subplot(515)
#plt.title('Rate of Climb ROC (m/s)')
#plt.plot(t,h,color='y',label='h')
plt.plot()