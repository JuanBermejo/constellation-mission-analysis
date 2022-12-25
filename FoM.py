import numpy as np
from TLE_OP_function import *
from read_celestrak import * 
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import CartesianRepresentation
from poliastro.util import time_range
from scipy.optimize import fsolve


# DATOS
Rt = 6378 # km
mu = 3.986e5 # km3/s2

# la GS está en Reyjkiavik (65.64737 N -20.24609 E) y se considera una elevación mínima de 15 deg
eps_GS = np.deg2rad(15) # rad
lat_GS = np.deg2rad(65.64737) # rad
long_GS = np.deg2rad(-20.24609) # rad

n_eval = 3000

constellation = (
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

# def visibility_status(constellation, lat_GS, long_GS, n_eval):

# Obtengo los parámetros orbitales para la última actualización del TLE del satélite
RAAN_deg, i_deg, e, a, AOP_deg, MA_deg = from_TLE_to_OrbParams(constellation)

# paso los datos de los parámetros orbitales a radianes para utilizar np.cos, np.sin, np.tan, etc.
RAAN =  np.deg2rad(RAAN_deg) # rad
i =  np.deg2rad(i_deg) # rad
AOP =  np.deg2rad(AOP_deg) # rad
MA =  np.deg2rad(MA_deg) # rad

# calculo el periodo de la órbita
T_sat = 2*np.pi*np.sqrt(a**3/mu) # s

# Considero que el satélite se encuentra en el perigeo en el momento en el que se ha actualizado el TLE
# Creo un vector de tiempos (en s) que va de 0 (paso por el perigeo) hasta el periodo orbital en n_eval tiempos
time_span = np.linspace(0, 30000, n_eval)

r_dot_cospsi = np.zeros(len(time_span))
w_rot_T = 2*np.pi/(24*3600) # rad/s

# Calculo el LST de la GS en el instante inicial
lat = 65.64737*u.deg
lon = -20.24609*u.deg
sat = list(load_gp_from_celestrak(name=constellation[0]))[0]
sat_dict = sat_to_dict(sat,constellation[0])
time = Time(sat_dict["EPOCH"], scale='utc')
LST_deg = time.sidereal_time('mean', longitude=lon)

for index, t in enumerate(time_span):
    # Hallo la anomalía excéntrica (rad) para cada instante de tiempo 
    def f(E):
        return E - e*np.sin(E) - np.sqrt(mu/a**3)*t
    E = fsolve(f, 0) # rad
    if E < 0:
        E = 2*np.pi + E
    # Hallo el radio (km) y la anomalía verdadera (rad)
    r = a*(1-e*np.cos(E)) # km
    TA = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2)) # rad
    if TA < 0:
        TA = 2*np.pi + TA
    # Hallo el tiempo sidereo local de la estación 
    LST_rad = LST_deg.value*np.pi/180 + w_rot_T*t
    # Ec(24) apartado 6.4 Elices 
    r_dot_cospsi[index] = r*( np.cos(lat_GS)*( np.cos(AOP+TA)*np.cos(LST_rad-RAAN) + np.cos(i)*np.sin(AOP+TA)*np.sin(LST_rad-RAAN) ) + np.sin(lat_GS)*np.sin(i)*np.sin(AOP+TA) )
    # r*( np.cos(lat_GS)*( np.cos(AOP+TA)*np.cos(LST_rad-RAAN) + np.cos(i)*np.sin(AOP+TA)*np.sin(LST_rad-RAAN) ) + np.sin(lat_GS)*np.sin(i)*np.sin(AOP+TA) )
    # print("r, r_dot_cospsi, E, TA, LST_rad, AOP, RAAN, i, lat_GS", r, r_dot_cospsi[index], E, TA, LST_rad, AOP, RAAN, i, lat_GS)

# Ec(21) apartado 6.4 Elices 
eps_sat = np.arcsin((r_dot_cospsi - Rt)/np.sqrt(r**2 + Rt**2 - 2*Rt*r_dot_cospsi))
# print(eps_sat*180/np.pi)
for j in range(len(eps_sat)):
    if eps_sat[j] < 0:
        eps_sat[j] = 2*np.pi + eps_sat[j]
# print(eps_sat*180/np.pi)

# Evalúo si la GS ve al satélite (True) o si no lo ve (False)
status = np.zeros(len(time_span))
for j in range(len(eps_sat)):
    if (eps_sat[j] >= eps_GS) and (eps_sat[j] <= np.pi- eps_GS):
        status[j] = 1 # true
    else:
        status[j] = 2 # false

print(time_span[np.where(status == np.amin(status))])
print(np.where(status == np.amin(status)) )
print(status[500])
print("t del 1er contacto =", time_span[523] - time_span[499])
print("t del 2o contacto =", time_span[1086] - time_span[1054])

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