# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 09:26:07 2018

@author: xx
"""

from numpy import*
from Cit_par import *



def rho00(hp):

    p=(1+Lambda*hp/Temp0)**(-g/(Lambda*R))*p0

    rho=(p/p0)**(1/gamma)*rho0
    return rho

def rho1(hp,Vc,Ttot): 
    T=Ttot-Vc**2/(2*cp)
    p=(1+Lambda*hp/Temp0)**(-g/(Lambda*R))*p0
    
    rho=p/(R*T)
    return rho