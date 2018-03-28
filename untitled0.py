# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 14:39:01 2018

@author: Max
"""

from Cit_par import *
from numpy import*
from control.matlab import*
import matplotlib.pyplot as plt
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

C1_symmetric=matrix([[-2*muc*c/(V0**2) ,0, 0 ,0],[0,(CZadot-2*muc)*c/V0,0,0],[0,0,-c/V0,0],[0,Cmadot*c/V0,0,-2*muc*KY2*c**2/V0**2]])
C2_symmetric=matrix([[CXu/V0,CXa,CZ0,CXq*c/V0],[CZu/V0,CZa,-CX0,(CZq+2*muc)*c/V0],[0,0,0,c/V0],[Cmu/V0,Cma,0,Cmq*c/V0]])
C3_symmetric=matrix([[CXde],[CZde],[0],[Cmde]])


A_symmetric=linalg.inv(-C1_symmetric)*C2_symmetric
B_symmetric=linalg.inv(-C1_symmetric)*C3_symmetric
C_symmetric=identity(4)
D_symmetric=zeros((4,1))

sys_symmetric=ss(A_symmetric,B_symmetric,C_symmetric,D_symmetric)

y,t,x = lsim(sys_symmetric,)