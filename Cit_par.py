# Citation 550 - Linear simulation
# xcg = 0.25 * c
from numpy import*

# Stationary flight condition

hp0    = 7000 * 0.3048 #change   	      # pressure altitude in the stationary flight condition [m]
V0     = 188.92* 0.5144444 #change               # true airspeed in the stationary flight condition [m/sec]
alpha0 = 0*pi/180 #change          # angle of attack in the stationary flight condition [rad]
th0    = 0*pi/180 #change          # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      = 60500/9.80665 #change          # mass [kg]

# aerodynamic properties
e      = 0.8772       # Oswald factor [ ]
CD0    = 0.0218        # Zero lift drag coefficient [ ]
CLa    = 4.64       # Slope of CL-alpha curve [ ]

# Longitudinal stability
Cma    = -0.666     # longitudinal stabilty [ ]
Cmde   = -1.345     # elevator effectiveness [ ]

# Aircraft geometry

S      = 30.00	          # wing area [m^2]
Sh     = 0.2 * S         # stabiliser area [m^2]
Sh_S   = Sh / S	          # [ ]
lh     = 0.71 * 5.968    # tail length [m]
c      = 2.0569	          # mean aerodynamic cord [m]
lh_c   = lh / c	          # [ ]
b      = 15.911	          # wing span [m]
bh     = 5.791	          # stabilser span [m]
A      = b ** 2 / S      # wing aspect ratio [ ]
Ah     = bh ** 2 / Sh    # stabilser aspect ratio [ ]
Vh_V   = 1	          # [ ]
ih     = -2 * pi / 180   # stabiliser angle of incidence [rad]

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

muc    = m / (rho * S * c)
mub    = m / (rho * S * b)
KX2    = 0.019
KZ2    = 0.042
KXZ    = 0.002
KY2    = 1.25 * 1.114

# Aerodynamic constants

Cmac   = 0                      # Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa                    # Wing normal force slope [ ]
CNha   = 2 * pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)            # Downwash gradient [ ]

# Lift and drag coefficient

CL = 2 * W / (rho * V0 ** 2 * S)              # Lift coefficient [ ]
CD = CD0 + (CLa * alpha0) ** 2 / (pi * A * e) # Drag coefficient [ ]

# Stabiblity derivatives

CX0    = W * sin(th0) / (0.5 * rho * V0 ** 2 * S)
CXu    = -0.08351
CXa    = 0
CXadot = +0.08330
CXq    = -0.28170
CXde   = -0.03728

CZ0    = -W * cos(th0) / (0.5 * rho * V0 ** 2 * S)
CZu    = -0.58115
CZa    = -5.74340
CZadot = -0.00350
CZq    = -5.66290
CZde   = -0.69612

Cmu    = +0.06990
Cmadot = +0.17800
Cmq    = -8.79415
Cm0    = 0.0297
CmTc   = -0.0064

CYb    = -0.7500
CYbdot =  0     
CYp    = -0.0304
CYr    = +0.8495
CYda   = -0.0400
CYdr   = +0.2300

Clb    = -0.10260
Clp    = -0.61707
Clr    = +0.23760
Clda   = -0.23088
Cldr   = +0.03440

Cnb    =  +0.12286
Cnbdot =   0     
Cnp    =  -0.0602
Cnr    =  -0.2073
Cnda   =  -0.0120
Cndr   =  -0.0939

