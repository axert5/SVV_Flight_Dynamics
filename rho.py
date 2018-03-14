# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 09:26:07 2018

@author: xx
"""

from math import*

p0=101325
T0=288.15
lamda=-0.0065
g=9.80665
R=287.05
gamma=1.4
rho0=1.225
cp=gamma*R/(gamma-1)


def rho(hp):

    p=(1+lamda*hp/T0)**(-g/(lamda*R))*p0

    rho=(p/p0)**(1/gamma)*rho0
    return rho,p

def rho2(hp,Vc): 
    T=Ttot-Vc**2/(2*cp)
    p=(1+lamda*hp/T0)**(-g/(lamda*R))*p0
    
    rho=p/(R*T)
    return rho,p


Vc=100
T=288.15
Ttot=T+Vc**2/(2*cp)
checklst=[]
for h in range(2000):
    check=rho(hp)[0]-rho2(hp,Vc)[0]
    checklst.append(check)
print (min(checklst),max(checklst))
    