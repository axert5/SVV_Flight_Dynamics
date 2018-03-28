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
from Flight_data import Short_period, Phugoid, Dutch_roll, Aperiodic_roll, Spiral
mass_start = 6702.89788

#Importing Flight data
Time_short, AoA_short, TAS_short, Pitch_short, Pitch_rate_short, Delta_e_short, Altitude_short, Fuel_left_short, Fuel_right_short = Short_period()
Time_phugoid, AoA_phugoid, TAS_phugoid, Pitch_phugoid, Pitch_rate_phugoid, Delta_e_phugoid, Altitude_phugoid, Fuel_left_phugoid, Fuel_right_phugoid = Phugoid()
Time_dutch, AoA_dutch, TAS_dutch, Pitch_dutch, Roll_dutch, Roll_rate_dutch, Yaw_rate_dutch, Delta_a_dutch, Delta_r_dutch, Altitude_dutch, Fuel_left_dutch, Fuel_right_dutch = Dutch_roll()
Time_aperiodic, AoA_aperiodic, TAS_aperiodic, Pitch_aperiodic, Roll_aperiodic, Roll_rate_aperiodic, Yaw_rate_aperiodic, Delta_a_aperiodic, Delta_r_aperiodic, Altitude_aperiodic, Fuel_left_aperiodic, Fuel_right_aperiodic = Aperiodic_roll()
Time_spiral, AoA_spiral, TAS_spiral, Pitch_spiral, Roll_spiral, Roll_rate_spiral, Yaw_rate_spiral, Delta_a_spiral, Delta_r_spiral, Altitude_spiral, Fuel_left_spiral, Fuel_right_spiral = Spiral()

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

#Initial dutch roll
changing_constant = changing_constants(Altitude_dutch[0], TAS_dutch[0], AoA_dutch[0], Pitch_dutch[0],(mass_start-Fuel_left_dutch[0]-Fuel_right_dutch[0]))
mucd = changing_constant[0]
mubd = changing_constant[1]
CLd  = changing_constant[2]
CDd  = changing_constant[3]
CX0d = changing_constant[4]
CZ0d = changing_constant[5]
V0d = TAS_dutch[0]

#Initial aperiodic
changing_constant = changing_constants(Altitude_aperiodic[0], TAS_aperiodic[0], AoA_aperiodic[0], Pitch_aperiodic[0],(mass_start-Fuel_left_aperiodic[0]-Fuel_right_aperiodic[0]))
muca = changing_constant[0]
muba = changing_constant[1]
CLa  = changing_constant[2]
CDa  = changing_constant[3]
CX0a = changing_constant[4]
CZ0a = changing_constant[5]
V0a = TAS_aperiodic[0]

#Initial spiral
changing_constant = changing_constants(Altitude_spiral[0], TAS_spiral[0], AoA_spiral[0], Pitch_spiral[0],(mass_start-Fuel_left_spiral[0]-Fuel_right_spiral[0]))
mucsp = changing_constant[0]
mubsp = changing_constant[1]
CLsp  = changing_constant[2]
CDsp  = changing_constant[3]
CX0sp = changing_constant[4]
CZ0sp = changing_constant[5]
V0sp = TAS_spiral[0]

#Correcting for initial deflections
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
Delta_a_dutch -= Delta_a_dutch[0]
Delta_r_dutch -= Delta_r_dutch[0]

TAS_spiral = TAS_spiral-TAS_spiral[0]
Roll_spiral = Roll_spiral-Roll_spiral[0]
Roll_rate_spiral -= Roll_rate_spiral[0] 
Delta_a_spiral -= Delta_a_spiral[0] - 0.0015
Delta_r_spiral -= Delta_r_spiral[0]

TAS_aperiodic = TAS_aperiodic-TAS_aperiodic[0]
Roll_aperiodic = Roll_aperiodic-Roll_aperiodic[0] 
Yaw_rate_aperiodic -= Yaw_rate_aperiodic[0]
Delta_a_aperiodic -= Delta_a_aperiodic[0]
Delta_r_aperiodic -= Delta_r_aperiodic[0]

#Optimisation of the model
mean_error = 1
plt.close('all')
while mean_error > 0.90:
#    CYb = random.uniform(-5,5)
#    CYp = random.uniform(-2,0)
#    CYr = random.uniform(-2,4)
#    Clb = random.uniform(-2,2)
#    Clp = random.uniform(-1,-0.4)
#    Clr = random.uniform(-1,1)
#    Cnb = random.uniform(0.0,0.18)
#    Cnp = random.uniform(-1,1)
#    Cnr = random.uniform(-0.3,-0.1)
    Cnr = -0.2073336615119732
    Clp = -0.6170747429900818
    Cnb = 0.12286105645574495
    
    #Setting the state space system dutch roll
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
    
    #Simulation
    u = np.transpose([Delta_a_dutch,-Delta_r_dutch])
    yd = lsim(sys_asymmetric, u, Time_dutch)
    
    #Error calculation
    Error = sqrt((1/len(Roll_dutch))*sum((yd[0][:,1]-Roll_dutch)**2))
    Relative_error_roll = Error/(max(Roll_dutch)-min(Roll_dutch))
    
    Error = sqrt((1/len(Roll_dutch))*sum((yd[0][:,2]-Roll_rate_dutch)**2))
    Relative_error_roll_rate = Error/(max(Roll_rate_dutch)-min(Roll_rate_dutch))
    
    Error = sqrt((1/len(Roll_dutch))*sum((yd[0][:,3]-Yaw_rate_dutch)**2))
    Relative_error_yaw_rate = Error/(max(Yaw_rate_dutch)-min(Yaw_rate_dutch))
    
    #Total average error
    Total_error =  Relative_error_roll + Relative_error_roll_rate + Relative_error_yaw_rate
    mean_errord = Total_error/3
    
    #Setting the state space system aperiodic roll
    C1_asymmetric = matrix([[(CYbdot-2*muba)*b/V0a,0,0,0],
                             [0,-b/(2*V0a),0,0], 
                             [0,0,(-4*muba*KX2*b**2)/(2*V0a**2),(4*muba*KXZ*b**2)/(2*V0a**2)],
                             [Cnbdot*b/(V0a),0,(4*muba*KXZ*b**2)/(2*V0a**2),(-4*muba*KZ2*b**2)/(2*V0a**2)]
                             ])
    
    C2_asymmetric = matrix([[CYb, CLa, CYp*b/(2*V0a),(CYr-4*muba)*b/(2*V0a)],
                             [0,0,b/(2*V0a),0],
                             [Clb,0,Clp*b/(2*V0a),Clr*b/(2*V0a)],
                             [Cnb,0,Cnp*b/(2*V0a),Cnr*b/(2*V0a)]
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
    
    #Simulation
    u = np.transpose([Delta_a_aperiodic,-Delta_r_aperiodic])
    ya = lsim(sys_asymmetric, u, Time_aperiodic)
    
    #Error Calculation aperiodic roll
    Error = sqrt((1/len(Roll_aperiodic))*sum((ya[0][:,1]-Roll_aperiodic)**2))
    Relative_error_roll = Error/(max(Roll_aperiodic)-min(Roll_aperiodic))
    
    Error = sqrt((1/len(Roll_aperiodic))*sum((ya[0][:,2]-Roll_rate_aperiodic)**2))
    Relative_error_roll_rate = Error/(max(Roll_rate_aperiodic)-min(Roll_rate_aperiodic))
    
    Error = sqrt((1/len(Roll_aperiodic))*sum((ya[0][:,3]-Yaw_rate_aperiodic)**2))
    Relative_error_yaw_rate = Error/(max(Yaw_rate_aperiodic)-min(Yaw_rate_aperiodic))
    
    #Total average error
    Total_error =  Relative_error_roll + Relative_error_roll_rate + Relative_error_yaw_rate
    mean_errora = Total_error/3
    
    #Setting the state space system spiral
    C1_asymmetric = matrix([[(CYbdot-2*mubsp)*b/V0sp,0,0,0],
                             [0,-b/(2*V0sp),0,0], 
                             [0,0,(-4*mubsp*KX2*b**2)/(2*V0sp**2),(4*mubsp*KXZ*b**2)/(2*V0sp**2)],
                             [Cnbdot*b/(V0sp),0,(4*mubsp*KXZ*b**2)/(2*V0sp**2),(-4*mubsp*KZ2*b**2)/(2*V0sp**2)]
                             ])
    
    C2_asymmetric = matrix([[CYb, CLsp, CYp*b/(2*V0sp),(CYr-4*mubsp)*b/(2*V0sp)],
                             [0,0,b/(2*V0sp),0],
                             [Clb,0,Clp*b/(2*V0sp),Clr*b/(2*V0sp)],
                             [Cnb,0,Cnp*b/(2*V0sp),Cnr*b/(2*V0sp)]
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
    
    #Simulation
    u = np.transpose([Delta_a_spiral, -Delta_r_spiral])
    ysp = lsim(sys_asymmetric, u, Time_spiral)
    
    #Error calculation spiral
    Error = sqrt((1/len(Roll_spiral))*sum((ysp[0][:,1]-Roll_spiral)**2))
    Relative_error_roll = Error/(max(Roll_spiral)-min(Roll_spiral))
    
    Error = sqrt((1/len(Roll_spiral))*sum((ysp[0][:,2]-Roll_rate_spiral)**2))
    Relative_error_roll_rate = Error/(max(Roll_rate_spiral)-min(Roll_rate_spiral))
    
    Error = sqrt((1/len(Roll_spiral))*sum((ysp[0][:,3]-Yaw_rate_spiral)**2))
    Relative_error_yaw_rate = Error/(max(Yaw_rate_spiral)-min(Yaw_rate_spiral))
    
    #Total average error
    Total_error =  Relative_error_roll + Relative_error_roll_rate + Relative_error_yaw_rate
    mean_errorsp = Total_error/3
    
    #Mean error of all aperiodic motions
    mean_error = (mean_errord + mean_errora + mean_errorsp)/3

print (mean_error)

#Plotting model
plt.figure('Dutch roll')
plt.subplot(221)
plt.plot(Time_dutch, Roll_dutch, label = 'Flight Data')
plt.plot(Time_dutch,yd[0][:,1], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('Roll angle [rad]')
plt.legend()


plt.subplot(222)
plt.plot(Time_dutch, Roll_rate_dutch, label = 'Flight Data')
plt.plot(Time_dutch,yd[0][:,2], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('Roll rate [rad/s]')
plt.legend()

plt.subplot(223)
plt.plot(Time_dutch, Yaw_rate_dutch, label = 'Flight Data')
plt.plot(Time_dutch,yd[0][:,3], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('Yaw rate [rad/s]')
plt.legend()
plt.show()

plt.subplot(224)
plt.plot(Time_dutch,Delta_a_dutch, label = 'Aileron')
plt.plot(Time_dutch,Delta_r_dutch, label = 'Rudder')
plt.xlabel('Time [s]')
plt.ylabel('Deflection [rad]')
plt.legend()
plt.show()


plt.figure('Aperiodic Roll')
plt.subplot(221)
plt.plot(Time_aperiodic, Roll_aperiodic, label = 'Flight Data')
plt.plot(Time_aperiodic,ya[0][:,1], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('Roll angle [rad]')
plt.legend()

plt.subplot(222)
plt.plot(Time_aperiodic, Roll_rate_aperiodic, label = 'Flight Data')
plt.plot(Time_aperiodic,ya[0][:,2], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('Roll rate [rad/s]')
plt.legend()

plt.subplot(223)
plt.plot(Time_aperiodic, Yaw_rate_aperiodic, label = 'Flight Data')
plt.plot(Time_aperiodic,ya[0][:,3], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('Yaw rate [rad/s]')
plt.legend()
plt.show()

plt.subplot(224)
plt.plot(Time_aperiodic,Delta_a_aperiodic, label = 'Aileron')
plt.plot(Time_aperiodic,Delta_r_aperiodic, label = 'Rudder')
plt.xlabel('Time [s]')
plt.ylabel('Deflection [rad]')
plt.legend()
plt.show()

plt.figure('Spiral')
plt.subplot(221)
plt.plot(Time_spiral, Roll_spiral, label = 'Flight Data')
plt.plot(Time_spiral,ysp[0][:,1], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('Roll angle [rad]')
plt.legend()

plt.subplot(222)
plt.plot(Time_spiral, Roll_rate_spiral, label = 'Flight Data')
plt.plot(Time_spiral,ysp[0][:,2], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('Roll rate [rad/s]')
plt.legend()

plt.subplot(223)
plt.plot(Time_spiral, Yaw_rate_spiral, label = 'Flight Data')
plt.plot(Time_spiral,ysp[0][:,3], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('Yaw rate [rad/s]')
plt.legend()
plt.show()

plt.subplot(224)
plt.plot(Time_spiral,Delta_a_spiral, label = 'Aileron')
plt.plot(Time_spiral,Delta_r_spiral, label = 'Rudder')
plt.xlabel('Time [s]')
plt.ylabel('Deflection [rad]')
plt.legend()
plt.show()

plt.figure('Symmetric inputs')
plt.subplot(121)
plt.plot(Time_short, Delta_e_short, label = 'Short Period')
plt.xlabel('Time [s]')
plt.ylabel('Deflection elevator [rad]')
plt.legend()

plt.subplot(122)
plt.plot(Time_phugoid, Delta_e_phugoid, label = 'Phugoid')
plt.xlabel('Time [s]')
plt.ylabel('Deflection elevator [rad]')
plt.legend()
plt.show()

mean_error = 1
while mean_error>0.90:
    Cmq = -6.883763192281982
    CZa = -5.832135330383927
    Cmadot = 0.12562411402000406
    CXu = -0.08398281071513852
    CXa = -0.03791622710921949
    CXq = -0.32156918147980135
    CZu = -0.4952377795039064
    CZadot = -0.007044212992164751
    CZq = -5.260621517613082
    Cmu = 0.07180053938491007
    
    #Cmq = random.uniform(0,-10)
    #Cmq = -4.753937622784239
    #CZa = random.uniform(-6,-3)
    #Cmadot = random.uniform(-0.5,0.5)
    #CXu = random.uniform(-0.5,0)
    #CXu = -0.11359455010228864
    #CXa = random.uniform(-1,0)
    #CXq = random.uniform(-0.4,-0.2)
    #CXq = -0.32977949703157144
    #CZu = random.uniform(-0.6,-0.2)
    #CZu = -0.4933392605640951
    #CZadot = random.uniform(-0.01,0)
    #CZq = random.uniform(-6,-4)
    #CZq = -5.628213819440996
    #Cmu = random.uniform(-0.1,0.1)
    #Cmu = 0.07653270911602936
    
    ##Setting the state space system short period
    C1_symmetric=matrix([[-2*mucs*c/(V0s**2) ,0, 0 ,0],[0,(CZadot-2*mucs)*c/V0s,0,0],[0,0,-c/V0s,0],[0,Cmadot*c/V0s,0,-2*mucs*KY2*c**2/V0s**2]])
    C2_symmetric=matrix([[CXu/V0s,CXa,CZ0s,CXq*c/V0s],[CZu/V0s,CZa,-CX0s,(CZq+2*mucs)*c/V0s],[0,0,0,c/V0s],[Cmu/V0s,Cma,0,Cmq*c/V0s]])
    C3_symmetric=matrix([[CXde],[CZde],[0],[Cmde]])
    
    
    A_symmetric=linalg.inv(-C1_symmetric)*C2_symmetric
    B_symmetric=linalg.inv(-C1_symmetric)*C3_symmetric
    C_symmetric=identity(4)
    D_symmetric=zeros((4,1))
    
    sys_symmetric=ss(A_symmetric,B_symmetric,C_symmetric,D_symmetric)
    
    #Simulation
    ys = lsim(sys_symmetric, Delta_e_short, Time_short)
    
    #Error calculation short period
    Error = sqrt((1/len(Pitch_rate_short))*sum((ys[0][:,0]-TAS_short)**2))
    Relative_error_TASs = Error/(max(TAS_short)-min(TAS_short))
    
    Error = sqrt((1/len(Pitch_rate_short))*sum((ys[0][:,1]-AoA_short)**2))
    Relative_error_AoAs = Error/(max(AoA_short)-min(AoA_short))
    
    Error = sqrt((1/len(Pitch_rate_short))*sum((ys[0][:,2]-Pitch_short)**2))
    Relative_error_Pitchs = Error/(max(Pitch_short)-min(Pitch_short))
    
    Error = sqrt((1/len(Pitch_rate_short))*sum((ys[0][:,3]-Pitch_rate_short)**2))
    Relative_error_Pitch_rates = Error/(max(Pitch_rate_short)-min(Pitch_rate_short))
    
    #Total average error
    Total_errors = Relative_error_TASs + Relative_error_AoAs + Relative_error_Pitchs + Relative_error_Pitch_rates
    mean_errors = Total_errors/4
    
    #Setting the state space system phugoid
    C1_symmetric=matrix([[-2*mucp*c/(V0p**2) ,0, 0 ,0],[0,(CZadot-2*mucp)*c/V0p,0,0],[0,0,-c/V0p,0],[0,Cmadot*c/V0p,0,-2*mucp*KY2*c**2/V0p**2]])
    C2_symmetric=matrix([[CXu/V0p,CXa,CZ0p,CXq*c/V0p],[CZu/V0p,CZa,-CX0p,(CZq+2*mucp)*c/V0p],[0,0,0,c/V0p],[Cmu/V0p,Cma,0,Cmq*c/V0p]])
    C3_symmetric=matrix([[CXde],[CZde],[0],[Cmde]])
    
    
    A_symmetric=linalg.inv(-C1_symmetric)*C2_symmetric
    B_symmetric=linalg.inv(-C1_symmetric)*C3_symmetric
    C_symmetric=identity(4)
    D_symmetric=zeros((4,1))
    
    Period_A_sym_dimless=2*pi/imag(array(linalg.eig(A_symmetric)[0]))
    print ('Period of Phugoid:',Period_A_sym_dimless[2])
    
    T12_A_sym_dimless=log(0.5)/real(array(linalg.eig(A_symmetric)[0]))
    print ('T1/2 of Phugoid:',T12_A_sym_dimless[2])
    
    sys_symmetric=ss(A_symmetric,B_symmetric,C_symmetric,D_symmetric)
    
    #Simulation
    yp = lsim(sys_symmetric, Delta_e_phugoid, Time_phugoid)
    
    #Error calculation phugoid
    Error = sqrt((1/len(Pitch_rate_phugoid))*sum((yp[0][:,0]-TAS_phugoid)**2))
    Relative_error_TASp = Error/(max(TAS_phugoid)-min(TAS_phugoid))
    
    Error = sqrt((1/len(Pitch_rate_phugoid))*sum((yp[0][:,1]-AoA_phugoid)**2))
    Relative_error_AoAp = Error/(max(AoA_phugoid)-min(AoA_phugoid))
    
    Error = sqrt((1/len(Pitch_rate_phugoid))*sum((yp[0][:,2]-Pitch_phugoid)**2))
    Relative_error_Pitchp = Error/(max(Pitch_phugoid)-min(Pitch_phugoid))
    
    Error = sqrt((1/len(Pitch_rate_phugoid))*sum((yp[0][:,3]-Pitch_rate_phugoid)**2))
    Relative_error_Pitch_ratep = Error/(max(Pitch_rate_phugoid)-min(Pitch_rate_phugoid))
    
    #Total average error
    Total_errorp = Relative_error_TASp + Relative_error_Pitchp + Relative_error_Pitch_ratep
    mean_errorp = Total_errorp/3
    
    #Mean error for all symmetric motions
    mean_error = (mean_errors+mean_errorp)/2 
    
print ('Mean error =', mean_error)

#Plotting graphs
plt.figure()
plt.plot(Time_phugoid,Delta_e_phugoid)
plt.show()
plt.figure('Short Period')
plt.subplot(221)
plt.plot(Time_short, TAS_short, label = 'Flight Data')
plt.plot(Time_short,ys[0][:,0], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('True Airspeed [m/s]')
plt.legend()


plt.subplot(222)
plt.plot(Time_short, AoA_short, label = 'Flight Data')
plt.plot(Time_short,ys[0][:,1], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('Angle of Attack [rad]')
plt.legend()

plt.subplot(223)
plt.plot(Time_short, Pitch_short, label = 'Flight Data')
plt.plot(Time_short,ys[0][:,2], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('Pitch angle [rad]')
plt.legend()

plt.subplot(224)
plt.plot(Time_short, Pitch_rate_short, label = 'Flight Data')
plt.plot(Time_short,ys[0][:,3], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('Pitch rate [rad/s]')
plt.legend()
plt.show()
    
plt.figure('Phugoid')
plt.subplot(221)
plt.plot(Time_phugoid, TAS_phugoid, label = 'Flight Data')
plt.plot(Time_phugoid,yp[0][:,0], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('True Airspeed [m/s]')
plt.legend()

plt.subplot(222)
plt.plot(Time_phugoid, AoA_phugoid, label = 'Flight Data')
plt.plot(Time_phugoid,yp[0][:,1], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('Angle of Attack [rad]')
plt.legend()

plt.subplot(223)
plt.plot(Time_phugoid, Pitch_phugoid, label = 'Flight Data')
plt.plot(Time_phugoid,yp[0][:,2], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('Pitch angle [rad]')
plt.legend()

plt.subplot(224)
plt.plot(Time_phugoid, Pitch_rate_phugoid, label = 'Flight Data')
plt.plot(Time_phugoid,yp[0][:,3], label = 'Simulation Data')
plt.xlabel('Time [s]')
plt.ylabel('Pitch rate [rad/s]')
plt.legend()
plt.show()
