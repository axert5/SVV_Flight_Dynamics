# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 09:40:12 2018

@author: Maria
"""

nr        = 6     #number of tests - should be the same for Excel and Python

excel_p    = [78178.80, 78178.80, 78178.80, 78178.80, 78208.49, 78535.73 ]
python_p   = [78178.23, 78178.23, 78178.23, 78178.20, 78207.92, 78535.17 ]
error_p    = []

excel_Vc   = [127.07, 112.66, 97.23, 81.80, 71.51, 59.16]
python_Vc  = [127.07, 112.66, 97.23, 81.80, 71.51, 59.16]
error_Vc   = []

excel_m    = [0.42, 0.38, 0.32, 0.27, 0.24, 0.20]
python_m   = [0.42, 0.38, 0.32, 0.27, 0.24, 0.20]
error_m    = []

excel_T    = [268.54, 268.09, 268.50, 268.94, 269.08, 268.56]
python_T   = [270.76, 268.09, 268.50, 268.94, 269.08, 268.50]
error_T    = []

excel_rho  = [1.01, 1.02, 1.01, 1.01, 1.01, 1.02]
python_rho = [1.00, 1.01, 1.01, 1.01, 1.01, 1.02]
error_rho  = []

excel_a    = [328.51, 328.24, 328.49, 328.76, 328.84, 328.53]
python_a   = [329.86, 328.23, 328.49, 328.75, 328.84, 328.52]
error_a    = []

excel_Vt   = [138.97, 123.24, 106.54, 89.78, 78.53, 64.81]
python_Vt  = [139.54, 123.24, 106.54, 89.78, 78.53, 64.80]
error_Vt   = []



for i in range (nr):
    
    if excel_p[i]<python_p[i]:
        error_p.append(100-(excel_p[i]/python_p[i]*100))
        
    else:
        error_p.append(100-(python_p[i]/excel_p[i]*100))
    
for i in range (nr):
    
    if excel_Vc[i]<python_Vc[i]:
        error_Vc.append(100-(excel_Vc[i]/python_Vc[i]*100))
        
    else:
        error_Vc.append(100-(python_Vc[i]/excel_Vc[i]*100))    

for i in range (nr):
    
    if excel_m[i]<python_m[i]:
        error_m.append(100-(excel_m[i]/python_m[i]*100))
        
    else:
        error_m.append(100-(python_m[i]/excel_m[i]*100))

for i in range (nr):
    
    if excel_T[i]<python_T[i]:
        error_T.append(100-(excel_T[i]/python_T[i]*100))
        
    else:
        error_T.append(100-(python_T[i]/excel_T[i]*100)) 
        
for i in range (nr):
    
    if excel_rho[i]<python_rho[i]:
        error_rho.append(100-(excel_rho[i]/python_rho[i]*100))
        
    else:
        error_rho.append(100-(python_rho[i]/excel_rho[i]*100))        
        
for i in range (nr):
    
    if excel_a[i]<python_a[i]:
        error_a.append(100-(excel_a[i]/python_a[i]*100))
        
    else:
        error_a.append(100-(python_a[i]/excel_a[i]*100)) 
        
for i in range (nr):
    
    if excel_Vt[i]<python_Vt[i]:
        error_Vt.append(100-(excel_Vt[i]/python_Vt[i]*100))
        
    else:
        error_Vt.append(100-(python_Vt[i]/excel_Vt[i]*100))        
        
"""        
print ("Errors on pressure \n", error_p, "\n max error on pressure: ", max(error_p) )        

print ("Errors on calibrated speed \n", error_Vc, "\n max error on calibrated airspeed: ", max(error_Vc) ) 

print ("Errors on mach \n", error_m, "\n max error on mach: ", max(error_m) ) 

print ("Errors on Temperature \n", error_T, "\n max error on temperature", max(error_T) ) 

print ("Errors on density \n ", error_rho, "\n max error on density", max(error_rho) )         

print ("Errors on speed of sound \n ", error_a, "\n max error on speed of sound", max(error_a) )         

print ("Errors on true airspeed \n", error_Vt, "\n max error on true airspeed", max(error_Vt) ) 
"""


print ("\n max error on pressure: ", max(error_p) )        

print ("\n max error on calibrated airspeed: ", max(error_Vc) ) 

print ("\n max error on mach: ", max(error_m) ) 

print ("\n max error on temperature", max(error_T) ) 

print ("\n max error on density", max(error_rho) )         

print ("\n max error on speed of sound", max(error_a) )         

print ("\n max error on true airspeed", max(error_Vt) ) 