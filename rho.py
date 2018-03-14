# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 09:26:07 2018

@author: xx
"""

from numpy import*

p0=101325
T0=288.15
lamda=-0.0065
g=9.80665
R=287.05
gamma=1.4
rho0=1.225
cp=gamma*R/(gamma-1)


def rho0(hp):

    p=(1+lamda*hp/T0)**(-g/(lamda*R))*p0

    rho=(p/p0)**(1/gamma)*rho0
    return rho

def rho(hp,Vc,Ttot): 
    T=Ttot-Vc**2/(2*cp)
    p=(1+lamda*hp/T0)**(-g/(lamda*R))*p0
    
    rho=p/(R*T)
    return rho