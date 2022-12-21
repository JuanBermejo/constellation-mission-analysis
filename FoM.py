import numpy as np
from TLE_OP_function import *
from read_celestrak import * 
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import CartesianRepresentation
from poliastro.util import time_range
from scipy.optimize import fsolve


FOSSA_list = (
    "FOSSASAT-2E2", 
    # "FOSSASAT-2E1",
    # "FOSSASAT-2E4", 
    # "FOSSASAT-2E3", 
    # "FOSSASAT-2E6", 
    # "FOSSASAT-2E5", 
    # "FOSSASAT-2E11",	
    # "FOSSASAT-2E12",	
    # "FOSSASAT-2E13",	
    # "FOSSASAT-2E7", 
    # "FOSSASAT-2E8" 
)

# DATOS
Rt = 6378 # km
mu = 3.986e5 # km3/s2
# la GS está en Reyjkiavik (65.64737 N -20.24609 E) y se considera una elevación mínima de 15 deg
eps_GS = np.deg2rad(15) # rad
lat_GS = np.deg2rad(65.64737) # rad
long_GS = np.deg2rad(-20.24609) # rad

RAAN_deg, i_deg, e, a, AOP_deg, MA_deg = from_TLE_to_OrbParams(FOSSA_list)
# paso los datos de los parámetros orbitales a radianes para utilizar np.cos, np.sin, np.tan, etc.
RAAN =  np.deg2rad(RAAN_deg) # rad
i =  np.deg2rad(i_deg) # rad
AOP =  np.deg2rad(AOP_deg) # rad
MA =  np.deg2rad(MA_deg) # rad
# calculo el periodo de la órbita
T_sat = 2*np.pi*np.sqrt(a**3/mu) # s

time_span = np.linspace(0, T_sat, 6) # vector de tiempos que va de 0 (paso por el perigeo) hasta el periodo orbital en 4 tiempos
r_dot_cospsi = np.zeros(len(time_span))
E = np.zeros(len(time_span))
TA = np.zeros(len(time_span))
r = np.zeros(len(time_span))

w_rot_T = 2*np.pi/(24*3600) # rad/s
GST0 = 0 # TBD!!!! 
for index, t in enumerate(time_span):
    # Hallo la anomalía excéntrica 
    def f(E):
        return E - e*np.sin(E) - np.sqrt(mu/a**3)*t
    E[index] = fsolve(f, 0) # rad
    if E[index] < 0:
        E[index] = 2*np.pi + E[index]
    # Hallo el radio y la anomalía verdadera
    r[index] = a*(1-e*np.cos(E[index])) # km
    TA[index] = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E[index]/2)) # rad
    if TA[index] < 0:
        TA[index] = 2*np.pi + TA[index]
    LST = GST0 + long_GS + w_rot_T*t
    r_dot_cospsi[index] = r[index]*( np.cos(lat_GS)*( np.cos(AOP+TA[index])*np.cos(LST-RAAN)+np.cos(i)*np.sin(AOP+TA[index])*np.sin(LST-RAAN) ) + np.sin(lat_GS)*np.sin(i)*np.sin(AOP+TA[index]) )

eps_sat = np.arcsin((r_dot_cospsi - Rt)/np.sqrt(r**2 + Rt**2 - 2*Rt*r_dot_cospsi))
print(eps_sat*180/np.pi)
for j in range(len(eps_sat)):
    if eps_sat[j] < 0:
        eps_sat[j] = 2*np.pi+eps_sat[j]
print(eps_sat*180/np.pi)

# ANTIGUO
# -------
# # Hallo la anomalía excéntrica y la anomalía verdadera
# def f(E):
#     return E - e*np.sin(E) - MA
# E = fsolve(f, 0) # rad
# TA = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2)) # rad
# if TA < 0:
#     TA = 2*np.pi + TA

# # Hallo el radio y el radiovector del satélite
# r_H = a*(1-e*np.cos(E)) # km
# r_H_vec =r_H*np.array([
#     np.cos(RAAN)*np.cos(AOP+TA)-np.sin(RAAN)*np.cos(i)*np.sin(AOP+TA), 
#     np.sin(RAAN)*np.cos(AOP+TA)+np.cos(RAAN)*np.cos(i)*np.sin(AOP+TA),
#     np.sin(i)*np.sin(AOP+TA)
# ]) # km

# LST = 0  
# c = np.array([np.cos(lambda_GS)*np.cos(LST), np.cos(lambda_GS)*np.sin(LST), np.sin(LST)])

# rho = np.zeros(len(r_H))
# eta = np.zeros(len(r_H))
# eps_sat = np.zeros(len(r_H))

# for j in range(len(r_H)):
#     rho[j] = np.arcsin( Rt/r_H[j] ) #rad
#     # DUDA: ARCTAN O ARCTAN2 !!
#     eta[j] = np.arctan( (np.sin(rho[j])*np.sin(lambda_GS))/(1-np.sin(rho[j])*np.cos(lambda_GS)) ) # rad
#     eps_sat[j] = np.cos(np.sin(eta[j])/np.sin(rho[j]))
# print(eps_sat*180/np.pi)

## DUDA
# sat = list(load_gp_from_celestrak(name=FOSSA_list[0]))[0]
# sat_dict = sat_to_dict(sat,FOSSA_list[0])
# now = Time(sat_dict["EPOCH"], scale='utc')
# now.jd1, now.jd2

# times = time_range(now, end=now + (T_sat/3600 << u.h), periods=3)

# errors, rs, vs = sat.sgp4_array(times.jd1, times.jd2)
# assert (errors == 0).all()
# rvec_CartesianRepresentation = CartesianRepresentation(rs << u.km, xyz_axis=-1)
# r_sat = rvec_CartesianRepresentation.norm().value
# # print("r_sat = ", r_sat)
# TA_sat = np.arccos( (a*(1-e**2)/r_sat - 1)/e ) 
# # print(TA, TA_sat)
# # print("rs =", rs)
# E_sat = np.arccos((1-r_sat/a)/e)
# print(E, E_sat)