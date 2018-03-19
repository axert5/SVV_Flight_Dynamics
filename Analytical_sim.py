"""
Program that gives values for each eigenmotion

"""

from Cit_par import *
from EigFuncVer import *
import numpy as np
from numpy.linalg import eig


"""
print("Dimension-less for short period are:" , short_period(muc , KY2 , CZa , Cmadot, Cmq, Cma))
print()
print("Dimension-less for phugoid are:" , phugoid(muc, CZa, Cmq, Cma,CXu , Cmu, CXa, CZu , CZ0))
print()
print("Dimension-less for Dutch roll are:" , dutch_roll(mub , KZ2 , Cnr , CYb , Cnb))
print ("-----------------------------------------------")
print("Dimension-having eigenvalues for short period are:" , short_period_d(muc , KY2 , CZa , Cmadot, Cmq, Cma , V0 , c))
print()
print("Dimension-having eigenvalues for phugoid are:" , phugoid_d(muc, CZa, Cmq, Cma,CXu , Cmu, CXa, CZu , CZ0, V0 , c))
print()
print("Dimension-having eigenvalues for Dutch roll are:" , dutch_roll_d(mub , KZ2 , Cnr , CYb , Cnb , V0 , b))
"""

#Taking real parts of each dimensionless eigenvalue
Re = np.array([[np.real(short_period(muc , KY2 , CZa , Cmadot, Cmq, Cma)[0])] ,
                [np.real(phugoid(muc, CZa, Cmq, Cma,CXu , Cmu, CXa, CZu , CZ0)[0])],
                [np.real(dutch_roll(mub , KZ2 , Cnr , CYb , Cnb)[0])]])

#Imaginary parts of each dimensionless eigenvalue
Im = np.array([[np.imag(short_period(muc , KY2 , CZa , Cmadot, Cmq, Cma)[0])] ,
                [np.imag(phugoid(muc, CZa, Cmq, Cma,CXu , Cmu, CXa, CZu , CZ0)[0])],
                [np.imag(dutch_roll(mub , KZ2 , Cnr , CYb , Cnb)[0])]])

#Short period data
P_sp       = 2*np.pi*c/(V0 * Im[0,0])                         #Period
omega0_sp  = V0 / c * np.sqrt(short_period(muc , KY2 , CZa , Cmadot, Cmq, Cma)[4] / short_period(muc , KY2 , CZa , Cmadot, Cmq, Cma)[2])     #Angular frequency of oscillation
damp_sp    = short_period(muc , KY2 , CZa , Cmadot, Cmq, Cma)[3] / (2 * np.sqrt(short_period(muc , KY2 , CZa , Cmadot, Cmq, Cma)[2] * short_period(muc , KY2 , CZa , Cmadot, Cmq, Cma)[4]))  #Damping coefficient
omegan_sp  = omega0_sp * np.sqrt(1-damp_sp**2)             #Eigenfrequency 
T_half_sp  = -np.log(0.5)  / omega0_sp / damp_sp                         #Half amplitude time

#Phugoid data
P_ph       = 2*np.pi*c/(V0 * Im[1,0])                         #Period
omega0_ph  = V0 / c * np.sqrt(phugoid(muc, CZa, Cmq, Cma,CXu , Cmu, CXa, CZu , CZ0)[4] / phugoid(muc, CZa, Cmq, Cma,CXu , Cmu, CXa, CZu , CZ0)[2])     #Angular frequency of oscillation
damp_ph    = -phugoid_damp(muc , CXu , CZu , CZ0)[3] / (2 * np.sqrt(phugoid_damp(muc , CXu , CZu , CZ0)[2] * phugoid_damp(muc , CXu , CZu , CZ0)[4]))  #Damping coefficient
omegan_ph  = omega0_ph * np.sqrt(1-damp_ph**2)             #Eigenfrequency 
T_half_ph  = -np.log(0.5) /omega0_ph / damp_ph            #Half amplitude time

#Dutch roll data
P_dr       = 2*np.pi*b/(V0 * Im[2,0])                         #Period
omega0_dr  = V0 / b * np.sqrt(dutch_roll(mub , KZ2 , Cnr , CYb , Cnb)[4] / dutch_roll(mub , KZ2 , Cnr , CYb , Cnb)[2])     #Angular frequency of oscillation
damp_dr    = dutch_roll(mub , KZ2 , Cnr , CYb , Cnb)[3] / (2 * np.sqrt(dutch_roll(mub , KZ2 , Cnr , CYb , Cnb)[2] * dutch_roll(mub , KZ2 , Cnr , CYb , Cnb)[4]))  #Damping coefficient
omegan_dr  = omega0_dr * np.sqrt(1-damp_dr**2)             #Eigenfrequency
T_half_dr  = -np.log(0.5)/ omega0_dr / damp_dr             #Half amplitude time



print ("------------Short period-------------")
print ("Period:             " , P_sp)
print ("T1/2 :              " , T_half_sp)
print ("Angular frequency:  " , omega0_sp)
print ("Damp coefficient:   " , damp_sp)
print ("Natural frequency:  " , omegan_sp) 
print()
print ("--------------Phugoid----------------")
print ("Period:             " , P_ph)
print ("T1/2 :              " , T_half_ph)
print ("Angular frequency:  " , omega0_ph)
print ("Damp coefficient:   " , damp_ph)
print ("Natural frequency:  " , omegan_ph) 
print()
print ("-------------Dutch roll--------------")
print ("Period:             " , P_dr)
print ("T1/2 :              " , T_half_dr)
print ("Angular frequency:  " , omega0_dr)
print ("Damp coefficient:   " , damp_dr)
print ("Natural frequency:  " , omegan_dr) 


A = 16 * mub**3 * (KX2 * KZ2 - KXZ**2)
B = -4 * mub**2 * (2 * CYb * (KX2*KZ2-KXZ**2) + Cnr * KX2 + Clp * KZ2 + (Clr + Cnp) * KXZ)
C = 2 * mub * ((CYb*Cnr - CYr * Cnb) * KX2 + (CYb * Clp - Clb * CYp) * KZ2 + 
               ((CYb * Cnp - Cnb * CYp) + (CYb * Clr - Clb * CYr)) * KXZ +
               4*mub * Cnb * KX2 + 4 * mub * Clb * KXZ + 0.5 * (Clp * Cnr - Cnp * Clr))
D = -4 * mub * CL * (Clb * KZ2 + Cnb * KXZ) + 2*mub*(Clb * Cnp - Cnb * Clp) + 0.5 * CYb * (Clr * Cnp - Cnr * Clp) + 0.5 * CYp * (Clb * Cnr- Cnb * Clr) + 0.5 * CYr * (Clp * Cnb - Cnp * Clb)
E = CL * (Clb * Cnr - Cnb * Clr)

Re = np.real(np.roots([A , B , C , D , E])[1])
Im = np.imag(np.roots([A , B , C , D , E])[1])

print ( -Re / np.sqrt(Re**2 + Im**2) , np.sqrt(Re**2 + Im**2) * V0 / b)