# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 14:14:30 2018

@author: xx

Weight and Balance
PH-LAB
Flight 1, 09.03.2017
"""
from numpy  import*
import matplotlib.pyplot as bplt
from scipy import interpolate

#inputs:
#0 Total current fuel mass of the aircraft
#1 cg shift. 0= no shift// 1= shift

#outputs:
#0 total aircraft mass in [kg]
#1 current x_cg from forward end of MAC in [m]
#2 total moment around forward end of MAC in [kgm]
#3 total moment around datum in [kgm]
#4 total aircraft mass in [lbs]
#5 current x_cg from datum in [in]
#6 total moment around datum in [in-lbs]


def cg(Fuel,shift):
    """
    inputs:
    0 Total current fuel mass of the aircraft
    1 cg shift. 0= no shift// 1= shift

    outputs:
    0 total aircraft mass in [kg]
    1 current x_cg from forward end of MAC in [m]
    2 total moment around forward end of MAC in [kgm]
    3 total moment around datum in [kgm]
    4 total aircraft mass in [lbs]
    5 current x_cg from datum in [in]
    6 total moment around datum in [in-lbs]
    """
    
    kg2lbs=1/0.453592
    in2m=0.0254
    
    BEW=9165.0
    xcg_BEW=292.18
    M_BEW=2677847.5 #hardcoded as above values are rounded. Most exact data used
    
    Seat_1=82*kg2lbs
    Seat_2=92*kg2lbs
    Seat_3=60*kg2lbs
    Seat_4=77*kg2lbs
    Seat_5=0*kg2lbs
    Seat_6=67*kg2lbs
    Seat_7=95*kg2lbs
    Seat_8=84*kg2lbs
    Seat_9=69*kg2lbs
    Seat_10=60*kg2lbs
    
    xcg_Seat12=131
    xcg_Seat34=214
    xcg_Seat56=251
    xcg_Seat78=288
    xcg_Seat910=170
#---------for cg shift, seat 7 moved to location between pilots    
    if shift==1:
        xcg_Seat7=131
    else:
        xcg_Seat7=xcg_Seat78
    
    M_Seat1=Seat_1*xcg_Seat12
    M_Seat2=Seat_2*xcg_Seat12
    M_Seat3=Seat_3*xcg_Seat34
    M_Seat4=Seat_4*xcg_Seat34
    M_Seat5=Seat_5*xcg_Seat56
    M_Seat6=Seat_6*xcg_Seat56
    M_Seat7=Seat_7*xcg_Seat7
    M_Seat8=Seat_8*xcg_Seat78
    M_Seat9=Seat_9*xcg_Seat910
    M_Seat10=Seat_10*xcg_Seat910
    
    
    ZFM=BEW+Seat_1+Seat_2+Seat_3+Seat_4+Seat_5+Seat_6+Seat_7+Seat_8+Seat_9+Seat_10
    M_ZFM=M_BEW+M_Seat1+M_Seat2+M_Seat3+M_Seat4+M_Seat5+M_Seat6+M_Seat7+M_Seat8+M_Seat9+M_Seat10
    xcg_ZFM=M_ZFM/ZFM
    
    #--------fuel-------------------------------------
    
    
    fuelmass=append(linspace(100,4900,49),array([5008]))
    fuelmoment=array([0298.16,0591.18,0879.08,1165.42,1448.40,1732.53,2014.80,2298.84,2581.92,
                      2866.30,3150.18,3434.52,3718.52,4003.23,4287.76,4572.24,4856.56,5141.16,
                      5425.64,5709.90,5994.04,6278.47,6562.82,6846.96,7131.00,7415.33,7699.60,
                      7984.34,8269.06,8554.05,8839.04,9124.80,9410.62,9696.97,9983.40,10270.08,
                      10556.84,10843.87,11131.00,11418.20,11705.50,11993.31,12281.18,12569.04,
                      12856.86,13144.73,13432.48,13720.56,14008.46,14320.34])
    
#total ramp fuel was 4100 lbs   
    min = 100
    u = 0
    for i in fuelmass:
        x=i-Fuel
        if x<0 and abs(x)< min:
            u = i
            min = x
    position= where(fuelmass == u)[0][0]
    
    M_fuel=(Fuel-fuelmass[position])/(fuelmass[position+1]-fuelmass[position])*(fuelmoment[position+1]-fuelmoment[position])+fuelmoment[position]        
    
    Current_mass=Fuel+ZFM
    Current_mass_kg=Current_mass/kg2lbs
    Current_moment=M_fuel+M_ZFM
    Current_moment_kgm=Current_moment/kg2lbs*in2m
    Current_xcg_datum=Current_moment/Current_mass
    Current_xcg_LEMAC=(Current_xcg_datum-261.45)*in2m
    Current_moment_LEMAC=Current_xcg_LEMAC*Current_mass_kg

    return Current_mass_kg,Current_xcg_LEMAC,Current_moment_LEMAC,Current_moment_kgm,Current_mass,Current_moment,Current_xcg_datum


        