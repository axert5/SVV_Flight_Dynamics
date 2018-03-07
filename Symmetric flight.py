# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 11:11:44 2018

@author: xx
"""
from Cit_par import *
#State space system

#symmetric case

C1_symmetric= matrix([[-2*muc*c/(V0**2) ,0, 0 ,0],[0,(CZadot-2*muc)*c/V0,0,0],[0,0,-c/V0,0],[0,Cma*c/V0,0,-2*muc*KY2**2*c**2/V0**2]])
C2_symmetric=matrix([[CXu/V0,CXa,CZ0,CXq*c/V0],[CZu/V0,CZa]])
