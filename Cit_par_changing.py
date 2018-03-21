# Citation 550 - Linear simulation
# xcg = 0.25 * c
from numpy import*

# Stationary flight condition
def changing_constants(hp0,V0,alpha0,th0,m):
    from Cit_par import e,CD0,CLa
    
    """
    Inputs:
        
    hp0       	 # pressure altitude in the stationary flight condition [m]
    V0           # true airspeed in the stationary flight condition [m/sec]
    alpha0       # angle of attack in the stationary flight condition [rad]
    th0          # pitch angle in the stationary flight condition [rad]
    m            # mass [kg]
    
    
    Outputs: muc,mub,CL,CD,CX0,CZ0
    
    muc = changing_constants[0]
    mub = changing_constants[1]
    CL  = changing_constants[2]
    CD  = changing_constants[3]
    CX0 = changing_constants[4]
    CZ0 = changing_constants[5]
    
    """
    # Constant values concerning atmosphere and gravity
    
    rho0   = 1.2250          # air density at sea level [kg/m^3] 
    Lambda = -0.0065         # temperature gradient in ISA [K/m]
    Temp0  = 288.15          # temperature at sea level in ISA [K]
    R      = 287.05          # specific gas constant [m^2/sec^2K]
    g      = 9.81            # [m/sec^2] (gravity constant)
    
    p0=101325
    gamma=1.4
    cp=gamma*R/(gamma-1)
    
    # air density [kg/m^3]  
    rho    = rho0 * power( ((1+(Lambda * hp0 / Temp0))), (-((g / (Lambda*R)) + 1)))   
    W      = m * g            # [N]       (aircraft weight)
    
    # Constant values concerning aircraft inertia
    S      = 30.00
    b      = 15.911
    c      = 2.0569
    A      = b ** 2 / S
       
    muc    = m / (rho * S * c)
    mub    = m / (rho * S * b)
    
    # Lift and drag coefficient
    
    CL = 2 * W / (rho * V0 ** 2 * S)              # Lift coefficient [ ]
    CD = CD0 + (CLa * alpha0) ** 2 / (pi * A * e) # Drag coefficient [ ]
    
    # Stabiblity derivatives
    
    CX0    = W * sin(th0) / (0.5 * rho * V0 ** 2 * S)
    CZ0    = -W * cos(th0) / (0.5 * rho * V0 ** 2 * S)
    return muc,mub,CL,CD,CX0,CZ0
