# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 16:49:12 2018

@author: xx
"""
import numpy as np

text_file=open("flightdata.txt","r")

line=np.loadtxt("flightdata.txt",delimiter=',',skiprows=1)
text_file.close()

begin_manoeuver=3497.0
position= where(line[:,47] == begin_manoeuver)[0]
#print (position)
print(line[:,13][34880])
print (line[:,14][34880])

print (line[:,41][34885])
"""
entry=entry_value-entry_vane_AOA
timestep=0.1s


vane_AOA
elevator_dte
column_fe
lh_engine_FMF
rh_engine_FMF
lh_engine_itt
rh_engine_itt
lh_engine_OP
rh_engine_OP
lh_engine_fan_N1
lh_engine_turbine_N2
rh_engine_fan_N1
rh_engine_turbine_N2
lh_engine_FU
rh_engine_FU
delta_a
delta_e
delta_r
Gps_date
Gps_utcSec
Ahrs1_Roll
Ahrs1_Pitch
Fms1_trueHeading
Gps_lat
Gps_long
Ahrs1_bRollRate
Ahrs1_bPitchRate
Ahrs1_bYawRate
Ahrs1_bLongAcc
Ahrs1_bLatAcc
Ahrs1_bNormAcc
Ahrs1_aHdgAcc
Ahrs1_xHdgAcc
Ahrs1_VertAcc
Dadc1_sat
Dadc1_tat
Dadc1_alt
Dadc1_bcAlt
Dadc1_bcAltMb
Dadc1_mach
Dadc1_cas
Dadc1_tas
Dadc1_altRate
measurement_running
measurement_n_rdy
display_graph_state
display_active_screen
time"""