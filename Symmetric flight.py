# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 11:11:44 2018

@author: xx
"""
from Cit_par import *
from numpy import*
from control.matlab import*
import matplotlib.pyplot as plt

#State space system
#testmaria
#symmetric case
print ('asdasdasd')
C1_symmetric= matrix([[-2*muc*c/(V0**2) ,0, 0 ,0],[0,(CZadot-2*muc)*c/V0,0,0],[0,0,-c/V0,0],[0,Cma*c/V0,0,-2*muc*KY2*c**2/V0**2]])
C2_symmetric=matrix([[CXu/V0,CXa,CZ0,CXq*c/V0],[CZu/V0,CZa,-CX0,(CZq+2*muc)*c/V0],[0,0,0,c/V0],[Cmu/V0,Cma,0,Cmq*c/V0]])
C3_symmetric=matrix([[CXde],[CZde],[0],[Cmde]])

#------original state space system---------------------------------------
A_symmetric=linalg.inv(-C1_symmetric)*C2_symmetric
B_symmetric=linalg.inv(-C1_symmetric)*C3_symmetric
C_symmetric=identity(4)
D_symmetric=zeros((4,1))

sys_symmetric=ss(A_symmetric,B_symmetric,C_symmetric,D_symmetric)

#print ('Eigenvalues of A_symmetric:',linalg.eig(A_symmetric)[0] )
#print ('Eigenvectors of A_symmetric:',linalg.eig(A_symmetric)[1] )


#------------state space expanded for altitude approximation----------------
temp_a=A_symmetric
temp_b=matrix([[0],[0],[0],[0]])
temp_c=matrix([[0,-V0,V0,0,0]])    #Assumption: Speed stays constant:V=V0; u=0=const
a_symmetric=vstack((hstack((temp_a, temp_b)),temp_c))
temp_d=B_symmetric
temp_e=matrix([[0]])
b_symmetric=vstack(((temp_d),temp_e))
c_symmetric=matrix([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1],[0,-V0,V0,0,0]])
d_symmetric=zeros((6,1))

sys_extended=ss(a_symmetric,b_symmetric,c_symmetric,d_symmetric)


#---state space computation-------------------------------------------------------

#-------inputs-------------
t=arange(0,10,0.01)
u=[-0.01]*1000

#--which model is selected----------------------------------------------------

#y=lsim(sys_symmetric,u,t) #original system


y=lsim(sys_extended,u,t) #system including altitude

# y[0][:,0]: u
# y[0][:,1]: alpha
# y[0][:,2]: theta
# y[0][:,3]: q
# y[0][:,4]: h   - for u=0
# y[0][:,5]: ROC - for u=0
#y[1]:       t

#----------plotting-----------------------------------------------------------
plt.figure(1)

plt.subplot(711)
plt.title('Input: Elevator Deflection Angle deltael (degs)')
plt.plot(t,array(u)*180/pi,color='m',label='i')

plt.subplot(712)
plt.title('Velocity u (m/s)')
plt.plot(t,y[0][:,0],color='c',label='u')

plt.subplot(713)
plt.title('Angle of Attack alpha (degs)')
plt.plot(t,y[0][:,1]*180/pi, color='r', label='alpha')

plt.subplot(714)
plt.title('Fligth Path Angle theta (degs)')
plt.plot(t,y[0][:,2]*180/pi,color='b',label='theta')

plt.subplot(715)
plt.title('Pitch Rate q (degs/s)')
plt.plot(t,y[0][:,3]*180/pi,color='g',label='q')

plt.subplot(716)
plt.title('Altitude (m)')
plt.plot(t,y[0][:,4],color='y',label='h')

plt.subplot(717)
plt.title('Rate of Climb (m/s)')
plt.plot(t,y[0][:,5],color='k',label='h')
plt.show()