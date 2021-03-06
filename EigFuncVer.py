"""
Code for the flight dynamics part of the SV0V0

Made by Andrei Badea

"""

import numpy as np

def short_period_d(muc , KY2 , CZa , Cmadot, Cmq, Cma , V0 , c):
    """This function outputs the  approximate dimension-having system eigenValues for short period oscillation.
    
    Inputs:
        muc     - relatiV0e density for symmetric motions m/(ro * S * c)
        Ky2     - non dimensional radius of gyration around Y axis (Ky2 / b, k = sqrt(Iy / m))
        CZa     - dCZ / da ; CZ = Z / (0.5 * ro * V0^2 * S) ; Z is the component of the aerodynamic force in the Z axis direction
        Cmadot  - dCm / d(adot * c / V0) ; Change in Cm with respect to change in adot * c / V0
        Cmq     - dCm / d(q * c / V0); Change in Cm with respect to q * c / V0
        Cma     - dCm / da ; Change in Cm wrt angle of attack
        V0       - Reference V0elocity
        c     - mean aerodynamic chord
                
    Outputs:
        Eig1    - First eigenValue of the system 
        Eig2    - Second eigenValue of the system 
    """
    A = 4 * muc**2 * KY2**2
    B = -2 * muc * (KY2**2 * CZa + Cmadot + Cmq)
    C = CZa * Cmq - 2 * muc * Cma
    j = np.complex(0 , 1)
    
    Eig1 = (-B + j * np.sqrt(4*A*C - B**2)) / 2 / A
    Eig2 = (-B - j * np.sqrt(4*A*C - B**2)) / 2 / A
    
    
    
    return Eig1 * V0 / c , Eig2 * V0 / c


def phugoid_d(muc, CZa, Cmq, Cma,CXu , Cmu, CXa, CZu , CZ0, V0 , c):
    """This function outputs the approximate dimension-having system eigenValues for phugoid motion.
    
    Inputs:
        muc     - relatiV0e density for symmetric motions m/(ro * S * c)
        CZa     - dCz / da ; Cz = Z / (0.5 * ro * V0^2 * S) ; Z is the component of the aerodynamic force in the Z axis direction
        Cmq     - dCm / d(q * c / V0); Change in Cm with respect to q * c / V0
        Cma     - dCm / da , change in moment coefficient with respect to angle of attack
        CXu     - dCX / d(q * c / V0) ; CX = X / (0.5 * ro * V0^2 * S); X is the component of the aerodynamic force in X direction
        Cmu     - dCm / du ; the change in Cm with respect to the change in the speed in X-axis "u"
        CXa     - dCX / da ; the change in the force in X direction with respect to angle of attack
        CZu     - dCZ / du ; the change in the force in Z direction with respect to the speed in "X" direction
        CZ0     - CZ in steady flight
        V0       - Reference V0elocity
        c     - mean aerodynamic chord
    
    Outputs:
        Eig1    - First eigenValue of the system
        Eig2    - Second eigenValue of the system
    """
    
    A = 2 * muc * (CZa * Cmq - 2 * muc * Cma)
    B = 2 * muc * (CXu * Cma - Cmu * CXa) + Cmq * (CZu * CXa - CXu * CZa)
    C = CZ0 * (Cmu * CZa - CZu * Cma)
    j = np.complex(0 , 1)
    
    Eig1 = (-B + j * np.sqrt(4*A*C - B**2)) / 2 / A
    Eig2 = (-B - j * np.sqrt(4*A*C - B**2)) / 2 / A
    
    return Eig1 * V0 / c, Eig2 * V0 / c

def dutch_roll_d(mub , KZ2 , Cnr , CYb , Cnb , V0 , b):
    """This function outputs the approximate dimension-having system eigenValues for the Dutch roll motion.
    
    Inputs:
        mub      - relatiV0e density for asymmetric motions m / (ro * S * b), b is the wing span
        KZ2      - non dimensional radius of gyration around the Z axis (KZ22 / b)
        Cnr     - dCn / d (rb/2V0); Cn is the yawing moment coefficient N / (0.5 * ro * V0^2 * S * b)
        CYb     - dCY / dbeta ; change in the force coefficient in Y direction wrt to change in the angle of sideslip (beta)
        Cnb     - dCn / dbeta ; Change in the yawing moment coeffienct wrt the change in the angle of sideslip (beta)
        V0       - Reference V0elocity
        b       - Wing span
    
    Outputs:
        Eig1    - First eigenValue of the system
        Eig2    - Second eigenValue of the system
    """
    A = 8 * mub**2 * KZ2 ** 2
    B = -2 * mub * (Cnr + 2 * KZ2**2 * CYb)
    C = 4 * mub * Cnb + CYb * Cnr
    j = np.complex(0 , 1)
    
    Eig1 = (-B + j * np.sqrt(4*A*C - B**2)) / 2 / A
    Eig2 = (-B - j * np.sqrt(4*A*C - B**2)) / 2 / A
    
    return Eig1 * V0 / b, Eig2 * V0  / b

def short_period(muc , KY2 , CZa , Cmadot, Cmq, Cma):
    """This function outputs the  approximate dimensionless system eigenValues for short period oscillation.
    
    Inputs:
        muc     - relatiV0e density for symmetric motions m/(ro * S * c)
        Ky2     - non dimensional radius of gyration around Y axis (Ky2 / b, k = sqrt(Iy / m))
        CZa     - dCZ / da ; CZ = Z / (0.5 * ro * V0^2 * S) ; Z is the component of the aerodynamic force in the Z axis direction
        Cmadot  - dCm / d(adot * c / V0) ; Change in Cm with respect to change in adot * c / V0
        Cmq     - dCm / d(q * c / V0); Change in Cm with respect to q * c / V0
        Cma     - dCm / da ; Change in Cm wrt angle of attack
                
    Outputs:
        Eig1    - First eigenValue of the system 
        Eig2    - Second eigenValue of the system 
    """
    A = 4 * muc**2 * KY2
    B = -2 * muc * (KY2 * CZa + Cmadot + Cmq)
    C = CZa * Cmq - 2 * muc * Cma
    j = np.complex(0 , 1)
    
    Eig1 = (-B + j * np.sqrt(4*A*C - B**2)) / 2 / A
    Eig2 = (-B - j * np.sqrt(4*A*C - B**2)) / 2 / A


    
    return Eig1  , Eig2 , A , B , C


def phugoid(muc, CZa, Cmq, Cma,CXu , Cmu, CXa, CZu , CZ0):
    """This function outputs the approximate dimensionless  system eigenValues for phugoid motion.
    
    Inputs:
        muc     - relatiV0e density for symmetric motions m/(ro * S * c)
        CZa     - dCz / da ; Cz = Z / (0.5 * ro * V0^2 * S) ; Z is the component of the aerodynamic force in the Z axis direction
        Cmq     - dCm / d(q * c / V0); Change in Cm with respect to q * c / V0
        Cma     - dCm / da , change in moment coefficient with respect to angle of attack
        CXu     - dCX / d(q * c / V0) ; CX = X / (0.5 * ro * V0^2 * S); X is the component of the aerodynamic force in X direction
        Cmu     - dCm / du ; the change in Cm with respect to the change in the speed in X-axis "u"
        CXa     - dCX / da ; the change in the force in X direction with respect to angle of attack
        CZu     - dCZ / du ; the change in the force in Z direction with respect to the speed in "X" direction
        CZ0     - CZ in steady flight

    
    Outputs:
        Eig1    - First eigenValue of the system
        Eig2    - Second eigenValue of the system
    """
    
    A = 2 * muc * (CZa * Cmq - 2 * muc * Cma)
    B = 2 * muc * (CXu * Cma - Cmu * CXa) + Cmq * (CZu * CXa - CXu * CZa)
    C = CZ0 * (Cmu * CZa - CZu * Cma)
    
    j = np.complex(0 , 1)
    
    Eig1 = (-B + j * np.sqrt(4*A*C - B**2)) / 2 / A
    Eig2 = (-B - j * np.sqrt(4*A*C - B**2)) / 2 / A
    
    return Eig1 , Eig2 , A , B , C

def dutch_roll(mub , KZ2 , Cnr , CYb , Cnb):
    """This function outputs the approximate dimensionless system eigenValues for the Dutch roll motion.
    
    Inputs:
        mub      - relatiV0e density for asymmetric motions m / (ro * S * b), b is the wing span
        KZ2      - non dimensional radius of gyration around the Z axis (KZ22 / b)
        Cnr     - dCn / d (rb/2V0); Cn is the yawing moment coefficient N / (0.5 * ro * V0^2 * S * b)
        CYb     - dCY / dbeta ; change in the force coefficient in Y direction wrt to change in the angle of sideslip (beta)
        Cnb     - dCn / dbeta ; Change in the yawing moment coeffienct wrt the change in the angle of sideslip (beta)

    
    Outputs:
        Eig1    - First eigenValue of the system
        Eig2    - Second eigenValue of the system
    """
    A = 8 * mub**2 * KZ2
    B = -2 * mub * (Cnr + 2 * KZ2 * CYb)
    C = 4 * mub * Cnb + CYb * Cnr
    j = np.complex(0 , 1)
    
    Eig1 = -B + j * np.sqrt(4*A*C - B**2) / 2 / A
    Eig2 = -B - j * np.sqrt(4*A*C - B**2) / 2 / A 
    
    return Eig1 , Eig2 , A , B , C

def phugoid_damp(muc , CXu , CZu , CZ0):
    A = -4*muc**2
    B = 2 * muc * CXu
    C = -CZu * CZ0
    j = np.complex(0 , 1)
    Eig1 = (-B + j * np.sqrt(4*A*C - B**2)) / 2 / A
    Eig2 = (-B - j * np.sqrt(4*A*C - B**2)) / 2 / A
    
    return Eig1 , Eig2 , A , B , C

def spiral(CL, Clb , Cnr , Cnb , Clr , Clp , CYb, Cnp , mub):
    
    Eig = (2*CL * (Clb * Cnr - Cnb * Clr)) / (Clp * (CYb * Cnr + 4 * mub * Cnb) - Cnp * (CYb * Clr + 4 * mub * Clb))
    
    return Eig

def aperiodic_roll(Clp , mub , KX2):
    
    Eig = Clp/(4 * mub * KX2)
    
    return Eig
    
    