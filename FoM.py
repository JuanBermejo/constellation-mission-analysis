import numpy as np
from TLE_OP_function import *
from read_celestrak import * 
from astropy import units as u
from astropy.time import Time
from poliastro.util import time_range
from scipy.optimize import fsolve


# DATOS
Rt = 6378 # km
mu = 3.986e5 # km3/s2

# la GS está en Reyjkiavik (65.64737 N -20.24609 E) y se considera una elevación mínima de 15 deg
eps_GS = np.deg2rad(15) # rad
lat_GS = np.deg2rad(65.64737) # rad
long_GS = np.deg2rad(-20.24609) # rad

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

# Considero que el satélite se encuentra en el perigeo en el momento en el que se ha actualizado el TLE
# Creo un vector de tiempos (en s) que va de 0 (paso por el perigeo) hasta 24 h evaluando casa segundo
t_fin = 23*3600 # s
delta_t = 2 # s
time_span = np.linspace(0, t_fin, int(t_fin/delta_t))

r_dot_cospsi = np.zeros(len(time_span))
w_rot_T = 2*np.pi/(24*3600) # rad/s

# Calculo el LST de la GS en el instante inicial
lat = 65.64737*u.deg
lon = -20.24609*u.deg
sat = list(load_gp_from_celestrak(name=constellation[0]))[0]
sat_dict = sat_to_dict(sat,constellation[0])
time = Time(sat_dict["EPOCH"], scale='utc')
LST_deg = time.sidereal_time('mean', longitude=lon) # h, min, s

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
    LST_rad = LST_deg.value*15*np.pi/180 + w_rot_T*t
    # Ec(24) apartado 6.4 Elices 
    r_dot_cospsi[index] = r*( np.cos(lat_GS)*( np.cos(AOP+TA)*np.cos(LST_rad-RAAN) + np.cos(i)*np.sin(AOP+TA)*np.sin(LST_rad-RAAN) ) + np.sin(lat_GS)*np.sin(i)*np.sin(AOP+TA) )

# Ec(21) apartado 6.4 Elices 
eps_sat = np.arcsin((r_dot_cospsi - Rt)/np.sqrt(r**2 + Rt**2 - 2*Rt*r_dot_cospsi))
for j in range(len(eps_sat)):
    if eps_sat[j] < 0:
        eps_sat[j] = 2*np.pi + eps_sat[j]

# Evalúo si la GS ve al satélite (True) o si no lo ve (False)
status = np.zeros(len(time_span))
for j in range(len(eps_sat)):
    if (eps_sat[j] >= eps_GS) and (eps_sat[j] <= np.pi- eps_GS):
        status[j] = 1 # true
    else:
        status[j] = 2 # false

contact = np.where(status == 1)[0]
contact_locator = []
contact_locator.append(contact[0])
for j in range(2,len(contact)):
    if contact[j]-contact[j-1] > 1:
        contact_locator.append(contact[j-1])
        contact_locator.append(contact[j])
contact_locator.append(contact[-1])

for j in range(0,len(contact_locator),2):
    print("Start time", Time(sat_dict["EPOCH"], format='isot', scale='utc') + time_span[contact_locator[j]]*u.s)
    print("Stop time", Time(sat_dict["EPOCH"], format='isot', scale='utc') + time_span[contact_locator[j+1]]*u.s)
    print("Duration in s of the contact", time_span[contact_locator[j+1]]-time_span[contact_locator[j]])
    print("------------------------------------------------------------------------------------------------")
print("Number of contacts =", len(contact_locator)/2)