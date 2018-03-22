# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 14:39:01 2018

@author: Max
"""

from Cit_par import *
from Cit_par_changing import *
from numpy import*
from control.matlab import*
import matplotlib.pyplot as plt
import warnings
import matplotlib.cbook
import random
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
from Flight_data import Short_period, Phugoid, Dutch_roll, Aperiodic_roll
mass_start = 6702.89788

Time_short, AoA_short, TAS_short, Pitch_short, Pitch_rate_short, Delta_e_short, Altitude_short, Fuel_left_short, Fuel_right_short = Short_period()
Time_phugoid, AoA_phugoid, TAS_phugoid, Pitch_phugoid, Pitch_rate_phugoid, Delta_e_phugoid, Altitude_phugoid, Fuel_left_phugoid, Fuel_right_phugoid = Phugoid()
Time_dutch, TAS_dutch, Roll_dutch, Roll_rate_dutch, Yaw_rate_dutch, Delta_a_dutch, Delta_r_dutch, Altitude_dutch, Fuel_left_dutch, Fuel_right_dutch = Dutch_roll()
Time_aperiodic, TAS_aperiodic, Roll_aperiodic, Roll_rate_aperiodic, Yaw_rate_aperiodic, Delta_a_aperiodic, Delta_r_aperiodic, Altitude_aperiodic, Fuel_left_aperiodic, Fuel_right_aperiodic = Aperiodic_roll()
#Initial short period
changing_constant = changing_constants(Altitude_short[0], TAS_short[0], AoA_short[0], Pitch_short[0],(mass_start-Fuel_left_short[0]-Fuel_right_short[0]))
mucs = changing_constant[0]
mubs = changing_constant[1]
CLs  = changing_constant[2]
CDs  = changing_constant[3]
CX0s = changing_constant[4]
CZ0s = changing_constant[5]
V0s = TAS_short[0]

#Initial phugoid
changing_constant = changing_constants(Altitude_phugoid[0], TAS_phugoid[0], AoA_phugoid[0], Pitch_phugoid[0],(mass_start-Fuel_left_phugoid[0]-Fuel_right_phugoid[0]))
mucp = changing_constant[0]
mubp = changing_constant[1]
CLp  = changing_constant[2]
CDp  = changing_constant[3]
CX0p = changing_constant[4]
CZ0p = changing_constant[5]
V0p = TAS_phugoid[0]

#Initial aperiodic roll
changing_constant = changing_constants(Altitude_phugoid[0], TAS_phugoid[0], AoA_phugoid[0], Pitch_phugoid[0],(mass_start-Fuel_left_phugoid[0]-Fuel_right_phugoid[0]))
mucd = changing_constant[0]
mubd = changing_constant[1]
CLd  = changing_constant[2]
CDd  = changing_constant[3]
CX0d = changing_constant[4]
CZ0d = changing_constant[5]
V0d = TAS_dutch[0]

TAS_short = TAS_short-TAS_short[0]
AoA_short = AoA_short-AoA_short[0]
Pitch_short = Pitch_short-Pitch_short[0]
Pitch_rate_short = Pitch_rate_short-Pitch_rate_short[0]

TAS_phugoid = TAS_phugoid-TAS_phugoid[0]
AoA_phugoid = AoA_phugoid-AoA_phugoid[0]
Pitch_phugoid = Pitch_phugoid-Pitch_phugoid[0]
Pitch_rate_phugoid = Pitch_rate_phugoid-Pitch_rate_phugoid[0]
 
TAS_dutch = TAS_dutch-TAS_dutch[0]
Roll_dutch = Roll_dutch-Roll_dutch[0] 
#plt.figure()
#plt.plot(Time_phugoid, Delta_e_phugoid)
#plt.show()
mean_error = 1
plt.figure()
plt.plot(Time_dutch, -Delta_a_dutch)
plt.plot(Time_dutch, -Delta_r_dutch)
plt.show()

while mean_error > 0.35:
    #CYb = random.uniform(-5,5)
    #CYp = random.uniform(-2,0)
    #CYr = random.uniform(-2,4)
    #Clb = random.uniform(-2,2)
    #Clp = random.uniform(-2,0)
    #Clr = random.uniform(-1,1)
    #Cnb = random.uniform(0,5)
    #Cnp = random.uniform(-1,1)
    #Cnr = random.uniform(-5,5)
    Cybdot = random.uniform(-5,5)
    Cnbdot = random.uniform(-5,5)
    #Cnr = -0.15
    #CYb = 0.6931206389752145
    #Cnb = 0.12362989451958462
    C1_asymmetric = matrix([[(CYbdot-2*mubd)*b/V0d,0,0,0],
                             [0,-b/(2*V0d),0,0], 
                             [0,0,(-4*mubd*KX2*b**2)/(2*V0d**2),(4*mubd*KXZ*b**2)/(2*V0d**2)],
                             [Cnbdot*b/(V0d),0,(4*mubd*KXZ*b**2)/(2*V0d**2),(-4*mubd*KZ2*b**2)/(2*V0d**2)]
                             ])
    
    C2_asymmetric = matrix([[CYb, CLd, CYp*b/(2*V0d),(CYr-4*mubd)*b/(2*V0d)],
                             [0,0,b/(2*V0d),0],
                             [Clb,0,Clp*b/(2*V0d),Clr*b/(2*V0d)],
                             [Cnb,0,Cnp*b/(2*V0d),Cnr*b/(2*V0d)]
                             ])
    
    C3_asymmetric = matrix([[CYda,CYdr],
                            [0,0],
                            [Clda,Cldr],
                            [Cnda,Cndr]
                            ])
    
    
    A_asymmetric = -linalg.inv(C1_asymmetric)*C2_asymmetric
    B_asymmetric = -linalg.inv(C1_asymmetric)*C3_asymmetric
    C_asymmetric = identity(4)
    D_asymmetric = zeros((4,2))
    
    sys_asymmetric=ss(A_asymmetric,B_asymmetric,C_asymmetric,D_asymmetric)
    
    u = np.transpose([-Delta_a_dutch,-Delta_r_dutch])
    y = lsim(sys_asymmetric, u, Time_dutch)
    
    Error = sqrt((1/len(Roll_dutch))*sum((y[0][:,1]-Roll_dutch)**2))
    Relative_error_roll = Error/(max(Roll_dutch)-min(Roll_dutch))
    #print ('TAS =',Relative_error_TAS)
    
    Error = sqrt((1/len(Roll_dutch))*sum((y[0][:,2]-Roll_rate_dutch)**2))
    Relative_error_roll_rate = Error/(max(Roll_rate_dutch)-min(Roll_rate_dutch))
    #print ('AoA =',Relative_error_AoA)
    
    Error = sqrt((1/len(Roll_dutch))*sum((y[0][:,3]-Yaw_rate_dutch)**2))
    Relative_error_yaw_rate = Error/(max(Yaw_rate_dutch)-min(Yaw_rate_dutch))
    #print ('Pitch =',Relative_error_Pitch)
    
    
    #Total average error
    Total_error =  Relative_error_roll + Relative_error_roll_rate + Relative_error_yaw_rate
    mean_error = Total_error/3

print(Relative_error_roll, Relative_error_roll_rate, Relative_error_yaw_rate)
print (mean_error)
plt.figure('Dutch roll')
plt.subplot(221)
plt.plot(Time_dutch,y[0][:,0])


plt.subplot(222)
plt.plot(Time_dutch, Roll_dutch)
plt.plot(Time_dutch,y[0][:,1])


plt.subplot(223)
plt.plot(Time_dutch, Roll_rate_dutch)
plt.plot(Time_dutch,y[0][:,2])


plt.subplot(224)
plt.plot(Time_dutch, Yaw_rate_dutch)
plt.plot(Time_dutch,y[0][:,3])
plt.show()

mean_error = 1
while mean_error > 0.20:
    #CYb = random.uniform(-5,5)
    #CYp = random.uniform(-2,0)
    #CYr = random.uniform(-2,4)
    #Clb = random.uniform(-2,2)
    #Clp = random.uniform(-2,0)
    #Clr = random.uniform(-1,1)
    #Cnb = random.uniform(0,5)
    #Cnp = random.uniform(-1,1)
    #Cnr = random.uniform(-5,5)
    #Cnr = -0.15
    #CYb = 0.6931206389752145
    #Cnb = 0.12362989451958462
    C1_asymmetric = matrix([[(CYbdot-2*mubd)*b/V0d,0,0,0],
                             [0,-b/(2*V0d),0,0], 
                             [0,0,(-4*mubd*KX2*b**2)/(2*V0d**2),(4*mubd*KXZ*b**2)/(2*V0d**2)],
                             [Cnbdot*b/(V0d),0,(4*mubd*KXZ*b**2)/(2*V0d**2),(-4*mubd*KZ2*b**2)/(2*V0d**2)]
                             ])
    
    C2_asymmetric = matrix([[CYb, CLd, CYp*b/(2*V0d),(CYr-4*mubd)*b/(2*V0d)],
                             [0,0,b/(2*V0d),0],
                             [Clb,0,Clp*b/(2*V0d),Clr*b/(2*V0d)],
                             [Cnb,0,Cnp*b/(2*V0d),Cnr*b/(2*V0d)]
                             ])
    
    C3_asymmetric = matrix([[CYda,CYdr],
                            [0,0],
                            [Clda,Cldr],
                            [Cnda,Cndr]
                            ])
    
    
    A_asymmetric = -linalg.inv(C1_asymmetric)*C2_asymmetric
    B_asymmetric = -linalg.inv(C1_asymmetric)*C3_asymmetric
    C_asymmetric = identity(4)
    D_asymmetric = zeros((4,2))
    
    sys_asymmetric=ss(A_asymmetric,B_asymmetric,C_asymmetric,D_asymmetric)
    
    u = np.transpose([-Delta_a_dutch,-Delta_r_dutch])
    y = lsim(sys_asymmetric, u, Time_dutch)
    
    Error = sqrt((1/len(Roll_dutch))*sum((y[0][:,1]-Roll_dutch)**2))
    Relative_error_roll = Error/(max(Roll_dutch)-min(Roll_dutch))
    #print ('TAS =',Relative_error_TAS)
    
    Error = sqrt((1/len(Roll_dutch))*sum((y[0][:,2]-Roll_rate_dutch)**2))
    Relative_error_roll_rate = Error/(max(Roll_rate_dutch)-min(Roll_rate_dutch))
    #print ('AoA =',Relative_error_AoA)
    
    Error = sqrt((1/len(Roll_dutch))*sum((y[0][:,3]-Yaw_rate_dutch)**2))
    Relative_error_yaw_rate = Error/(max(Yaw_rate_dutch)-min(Yaw_rate_dutch))
    #print ('Pitch =',Relative_error_Pitch)
    
    
    #Total average error
    Total_error =  Relative_error_roll_rate + Relative_error_yaw_rate
    mean_error = Total_error/2

print(Relative_error_roll, Relative_error_roll_rate, Relative_error_yaw_rate)
print (mean_error)
plt.figure('Dutch roll 2')
plt.subplot(221)
plt.plot(Time_dutch,y[0][:,0])


plt.subplot(222)
plt.plot(Time_dutch, Roll_dutch)
plt.plot(Time_dutch,y[0][:,1])


plt.subplot(223)
plt.plot(Time_dutch, Roll_rate_dutch)
plt.plot(Time_dutch,y[0][:,2])


plt.subplot(224)
plt.plot(Time_dutch, Yaw_rate_dutch)
plt.plot(Time_dutch,y[0][:,3])
plt.show()
'''
mean_error = 1
while mean_error>0.90:
    
    Cmq = random.uniform(0,-10)
    #Cmq = -4.753937622784239
    CZa = random.uniform(-6,-3)
    Cmadot = random.uniform(-0.5,0.5)
    CXu = random.uniform(-0.5,0)
    #CXu = -0.11359455010228864
    CXa = random.uniform(-1,0)
    CXq = random.uniform(-0.4,-0.2)
    #CXq = -0.32977949703157144
    CZu = random.uniform(-0.6,-0.2)
    #CZu = -0.4933392605640951
    CZadot = random.uniform(-0.01,0)
    CZq = random.uniform(-6,-4)
    #CZq = -5.628213819440996
    Cmu = random.uniform(-0.1,0.1)
    #Cmu = 0.07653270911602936
    
    C1_symmetric=matrix([[-2*mucs*c/(V0s**2) ,0, 0 ,0],[0,(CZadot-2*mucs)*c/V0s,0,0],[0,0,-c/V0s,0],[0,Cmadot*c/V0s,0,-2*mucs*KY2*c**2/V0s**2]])
    C2_symmetric=matrix([[CXu/V0s,CXa,CZ0s,CXq*c/V0s],[CZu/V0s,CZa,-CX0s,(CZq+2*mucs)*c/V0s],[0,0,0,c/V0s],[Cmu/V0s,Cma,0,Cmq*c/V0s]])
    C3_symmetric=matrix([[CXde],[CZde],[0],[Cmde]])
    
    
    A_symmetric=linalg.inv(-C1_symmetric)*C2_symmetric
    B_symmetric=linalg.inv(-C1_symmetric)*C3_symmetric
    C_symmetric=identity(4)
    D_symmetric=zeros((4,1))
    
    sys_symmetric=ss(A_symmetric,B_symmetric,C_symmetric,D_symmetric)
    
    ys = lsim(sys_symmetric, Delta_e_short, Time_short)
    
    Error = sqrt((1/len(Pitch_rate_short))*sum((ys[0][:,0]-TAS_short)**2))
    Relative_error_TASs = Error/(max(TAS_short)-min(TAS_short))
    #print ('TAS =',Relative_error_TAS)
    
    Error = sqrt((1/len(Pitch_rate_short))*sum((ys[0][:,1]-AoA_short)**2))
    Relative_error_AoAs = Error/(max(AoA_short)-min(AoA_short))
    #print ('AoA =',Relative_error_AoA)
    
    Error = sqrt((1/len(Pitch_rate_short))*sum((ys[0][:,2]-Pitch_short)**2))
    Relative_error_Pitchs = Error/(max(Pitch_short)-min(Pitch_short))
    #print ('Pitch =',Relative_error_Pitch)
    
    Error = sqrt((1/len(Pitch_rate_short))*sum((ys[0][:,3]-Pitch_rate_short)**2))
    Relative_error_Pitch_rates = Error/(max(Pitch_rate_short)-min(Pitch_rate_short))
    #print ('Pitch rate =',Relative_error_Pitch_rate)
    
    #Total average error
    Total_errors = Relative_error_TASs + Relative_error_AoAs + Relative_error_Pitchs + Relative_error_Pitch_rates
    mean_errors = Total_errors/4
    
    #plt.figure()
    #plt.plot(Time_short, Delta_e_short)
    #plt.show()
    
    #print (Altitude_phugoid[0], TAS_phugoid[0], AoA_phugoid[0], Pitch_phugoid[0],(mass_start-Fuel_left_phugoid[0]-Fuel_right_phugoid[0]))
    #Initial Phugoid
    
    

    
    C1_symmetric=matrix([[-2*mucp*c/(V0p**2) ,0, 0 ,0],[0,(CZadot-2*mucp)*c/V0p,0,0],[0,0,-c/V0p,0],[0,Cmadot*c/V0p,0,-2*mucp*KY2*c**2/V0p**2]])
    C2_symmetric=matrix([[CXu/V0p,CXa,CZ0p,CXq*c/V0p],[CZu/V0p,CZa,-CX0p,(CZq+2*mucp)*c/V0p],[0,0,0,c/V0p],[Cmu/V0p,Cma,0,Cmq*c/V0p]])
    C3_symmetric=matrix([[CXde],[CZde],[0],[Cmde]])
    
    
    A_symmetric=linalg.inv(-C1_symmetric)*C2_symmetric
    B_symmetric=linalg.inv(-C1_symmetric)*C3_symmetric
    C_symmetric=identity(4)
    D_symmetric=zeros((4,1))
    
    sys_symmetric=ss(A_symmetric,B_symmetric,C_symmetric,D_symmetric)
    
    yp = lsim(sys_symmetric, Delta_e_phugoid, Time_phugoid)
    
    Error = sqrt((1/len(Pitch_rate_phugoid))*sum((yp[0][:,0]-TAS_phugoid)**2))
    Relative_error_TASp = Error/(max(TAS_phugoid)-min(TAS_phugoid))
    #print ('TAS =',Relative_error_TAS)
    
    Error = sqrt((1/len(Pitch_rate_phugoid))*sum((yp[0][:,1]-AoA_phugoid)**2))
    Relative_error_AoAp = Error/(max(AoA_phugoid)-min(AoA_phugoid))
    #print ('AoA =',Relative_error_AoA)
    
    Error = sqrt((1/len(Pitch_rate_phugoid))*sum((yp[0][:,2]-Pitch_phugoid)**2))
    Relative_error_Pitchp = Error/(max(Pitch_phugoid)-min(Pitch_phugoid))
    #print ('Pitch =',Relative_error_Pitch)
    
    Error = sqrt((1/len(Pitch_rate_phugoid))*sum((yp[0][:,3]-Pitch_rate_phugoid)**2))
    Relative_error_Pitch_ratep = Error/(max(Pitch_rate_phugoid)-min(Pitch_rate_phugoid))
    #print ('Pitch rate =',Relative_error_Pitch_rate)
    
    #Total average error
    Total_errorp = Relative_error_TASp + Relative_error_Pitchp + Relative_error_Pitch_ratep
    mean_errorp = Total_errorp/3
    
    mean_error = (mean_errors+mean_errorp)/2
    #mean_error = Relative_error_AoAp   
    
    
print (Relative_error_TASs, Relative_error_AoAs, Relative_error_Pitchs, Relative_error_Pitch_rates)
print (Relative_error_TASp, Relative_error_AoAp, Relative_error_Pitchp, Relative_error_Pitch_ratep)
print ('Short AoA =', Relative_error_AoAs)
print ('Phugoid AoA =', Relative_error_AoAp)
print ('Mean error =', mean_error)
plt.figure('Short Period')
plt.subplot(221)
plt.plot(Time_short, TAS_short)
plt.plot(Time_short,ys[0][:,0])


plt.subplot(222)
plt.plot(Time_short, AoA_short)
plt.plot(Time_short,ys[0][:,1])


plt.subplot(223)
plt.plot(Time_short, Pitch_short)
plt.plot(Time_short,ys[0][:,2])


plt.subplot(224)
plt.plot(Time_short, Pitch_rate_short)
plt.plot(Time_short,ys[0][:,3])
plt.show()
    
plt.figure('Phugoid')
plt.subplot(221)
plt.plot(Time_phugoid, TAS_phugoid)
plt.plot(Time_phugoid,yp[0][:,0])


plt.subplot(222)
plt.plot(Time_phugoid, AoA_phugoid)
plt.plot(Time_phugoid,yp[0][:,1])


plt.subplot(223)
plt.plot(Time_phugoid, Pitch_phugoid)
plt.plot(Time_phugoid,yp[0][:,2])


plt.subplot(224)
plt.plot(Time_phugoid, Pitch_rate_phugoid)
plt.plot(Time_phugoid,yp[0][:,3])
plt.show()
'''
##Relative Root Mean Square Error
#Error = sqrt((1/len(Pitch_rate_short))*sum((y[0][:,0]-TAS_short)**2))
#Relative_error_TAS = Error/(max(TAS_short)-min(TAS_short))
#print ('TAS =',Relative_error_TAS)
#
#Error = sqrt((1/len(Pitch_rate_short))*sum((y[0][:,1]-AoA_short)**2))
#Relative_error_AoA = Error/(max(AoA_short)-min(AoA_short))
#print ('AoA =',Relative_error_AoA)
#
#Error = sqrt((1/len(Pitch_rate_short))*sum((y[0][:,2]-Pitch_short)**2))
#Relative_error_Pitch = Error/(max(Pitch_short)-min(Pitch_short))
#print ('Pitch =',Relative_error_Pitch)
#
#Error = sqrt((1/len(Pitch_rate_short))*sum((y[0][:,3]-Pitch_rate_short)**2))
#Relative_error_Pitch_rate = Error/(max(Pitch_rate_short)-min(Pitch_rate_short))
#print ('Pitch rate =',Relative_error_Pitch_rate)
#
##Total average error
#Total_error = Relative_error_TAS + Relative_error_AoA + Relative_error_Pitch + Relative_error_Pitch_rate
#Mean_error = Total_error/4
#print ('Mean total error =', Mean_error)