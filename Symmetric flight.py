# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 11:11:44 2018

@author: xx
"""
from Cit_par_book import *
from Cit_par_changing import changing_constants
from WeightBalance import cg
from numpy import*
from control.matlab import*
import matplotlib.pyplot as plt
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)


#symmetric case
'''
#begin short period
hp0     =7000*0.3048
V0      =188.92*0.51444
alpha0  =5*pi/180
th0     =0*pi/180
fuel_used_LEngine=504.276941303232
fuel_used_REngine=520.152402665631
m       =cg(4100-(fuel_used_LEngine+fuel_used_REngine),0)[0]


changing_constants=changing_constants(hp0,V0,alpha0,th0,m)

muc = changing_constants[0]
mub = changing_constants[1]
CL  = changing_constants[2]
CD  = changing_constants[3]
CX0 = changing_constants[4]
CZ0 = changing_constants[5]'''
#------original state space system---------------------------------------

#dimension having
C1_symmetric=matrix([[-2*muc*c/(V0**2) ,0, 0 ,0],[0,(CZadot-2*muc)*c/V0,0,0],[0,0,-c/V0,0],[0,Cmadot*c/V0,0,-2*muc*KY2*c**2/V0**2]])
C2_symmetric=matrix([[CXu/V0,CXa,CZ0,CXq*c/V0],[CZu/V0,CZa,-CX0,(CZq+2*muc)*c/V0],[0,0,0,c/V0],[Cmu/V0,Cma,0,Cmq*c/V0]])
C3_symmetric=matrix([[CXde],[CZde],[0],[Cmde]])


A_symmetric=linalg.inv(-C1_symmetric)*C2_symmetric
B_symmetric=linalg.inv(-C1_symmetric)*C3_symmetric
C_symmetric=identity(4)
D_symmetric=zeros((4,1))

sys_symmetric=ss(A_symmetric,B_symmetric,C_symmetric,D_symmetric)

#dimensionless
C1_sym_dimless=matrix([[-2*muc*c/V0,0,0,0],[0,(CZadot-2*muc)*c/V0,0,0],[0,0,-c/V0,0],[0,Cmadot*c/V0,0,-2*muc*KY2*c/V0]])
C2_sym_dimless=matrix([[CXu,CXa,CZ0,CXq],[CZu,CZa,CX0,(CZq+2*muc)],[0,0,0,1],[Cmu,Cma,0,Cmq]])
C3_sym_dimless=matrix([[CXde],[CZde],[0],[Cmde]])

A_sym_dimless=linalg.inv(-C1_sym_dimless)*C2_sym_dimless
B_sym_dimless=linalg.inv(-C1_sym_dimless)*C3_sym_dimless
C_sym_dimless=identity(4)
C_sym_hybrid=matrix([[V0,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,V0/c]])
D_sym_dimless=zeros((4,1))

#true dimensionless output
sys_sym_dimless=ss(A_sym_dimless,B_sym_dimless,C_sym_dimless,D_sym_dimless)

#dimensionless computation but dimension having output
sys_sym_hybrid=ss(A_sym_dimless,B_sym_dimless,C_sym_hybrid,D_sym_dimless)

#------------state space expanded for altitude approximation----------------
temp_a=A_sym_dimless
temp_b=matrix([[0],[0],[0],[0]])
temp_c=matrix([[0,-V0,V0,0,0]])    #Assumption: Speed stays constant:V=V0; u=0=const
a_sym_dimless=vstack((hstack((temp_a, temp_b)),temp_c))
temp_d=B_sym_dimless
temp_e=matrix([[0]])
b_sym_dimless=vstack(((temp_d),temp_e))
c_sym_dimless=matrix([[V0,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,V0/c,0],[0,0,0,0,1],[0,-V0,V0,0,0]])
d_sym_dimless=zeros((6,1))

sys_extended=ss(a_sym_dimless,b_sym_dimless,c_sym_dimless,d_sym_dimless)


#---state space computation------------------------------------------------------
#-------------------------------------------------------------------------------

#-------inputs-------------
t=arange(0,15,0.01)

ude=[-0.012]*len(t) #input vector for elevator deflection

#--which model is selected----------------------------------------------------

#sys=sys_symmetric           #standard dimension having
#sys=sys_sym_hybrid          #dimless computation, dim having outputs
sys=sys_sym_dimless         #dimless outputs
#sys=sys_extended            #dimension having, extended for approx. ROC and altitude       


#--what Input is selected------------------------------------------------------



y=lsim(sys,ude,t)   #Using input vector
#y=impulse(sys,t)                #Impulse input
#y=step(sys,t)                   #Step input

# y[0][:,0]: u
# y[0][:,1]: alpha
# y[0][:,2]: theta
# y[0][:,3]: q
# y[0][:,4]: h   - for u=0
# y[0][:,5]: ROC - for u=0
#y[1]:       t

#-----Matrix analysis--------------------------------------------------------

#--------dimensionless system--------------------------------------------------
print()

print('Symmetric Flight:')

print()

eigenvalues_A_sym_dimless=linalg.eig(A_sym_dimless)[0]
print ('Eigenvalues of Short Period:',eigenvalues_A_sym_dimless[:-2] )
print ('Eigenvalues of Phugoid:',eigenvalues_A_sym_dimless[-2:] )

print()

T12_A_sym_dimless=log(0.5)/real(array(linalg.eig(A_sym_dimless)[0]))
print ('T1/2 of Short Period:',T12_A_sym_dimless[0])
print ('T1/2 of Phugoid:',T12_A_sym_dimless[2])

print()

Period_A_sym_dimless=2*pi/imag(array(linalg.eig(A_sym_dimless)[0]))
print ('Period of Short Period:',Period_A_sym_dimless[0])
print ('Period of Phugoid:',Period_A_sym_dimless[2])

print()

damping=damp(sys,doprint=False)[1]
print ('Damping of Short Period:',damping[0])
print ('Damping of Phugoid:',damping[2])

print()

natfreq=damp(sys,doprint=False)[0]*sqrt(1-damp(sys,doprint=False)[1]**2)
print ('Nat. Frequency of Short Period:',natfreq[0])
print ('Nat. Frequency of Phugoid:',natfreq[2])

print()

print ('Frequency of Short Period:',natfreq[0] / sqrt(1 - damping[0]**2))
print ('Frequency of Phugoid:',natfreq[2]/ sqrt(1 - damping[2]**2))

print()
#----------plotting-----------------------------------------------------------
plt.figure(1)

plt.subplot(511)
plt.title('Input: Elevator Deflection Angle deltael (rad)')
plt.plot(t,array(ude),color='m',label='i')

plt.subplot(512)
plt.title('Velocity u (m/s)')
plt.plot(t,y[0][:,0],color='c',label='u')

plt.subplot(513)
plt.title('Angle of Attack alpha (rad)')
plt.plot(t,y[0][:,1], color='r', label='alpha')

plt.subplot(514)
plt.title('Fligth Path Angle theta (rad)')
plt.plot(t,y[0][:,2],color='b',label='theta')

plt.subplot(515)
plt.title('Pitch Rate q (rad/s)')
plt.plot(t,y[0][:,3],color='g',label='q')
"""
plt.subplot(716)
plt.title('Altitude (m)')
#plt.plot(t,y[0][:,4],color='y',label='h')

plt.subplot(717)
plt.title('Rate of Climb (m/s)')
#plt.plot(t,y[0][:,5],color='k',label='h')"""

#plt.show()