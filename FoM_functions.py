import numpy as np
import pandas as pd
from astropy import units as u
import csv

from TLE_OP_function import from_TLE_to_OrbParams, E_and_TA_from_MA, time_in_sat_epoch, time_in_constellation_epoch
from scipy.optimize import fsolve

# DATOS
Rt = 6378 # km
mu = 3.986e5 # km3/s2

def elevacion_sat(time_span, lat_GS, long_GS, constellation):
    # Obtengo los parámetros orbitales para la última actualización del TLE del satélite
    RAAN_deg, i_deg, e, a, AOP_deg, MA_deg = from_TLE_to_OrbParams(constellation)    # dim (sat, )

    # paso los datos de los parámetros orbitales a radianes para utilizar np.cos, np.sin, np.tan, etc.
    RAAN =  np.deg2rad(RAAN_deg) # rad, dim (sat, )
    i =  np.deg2rad(i_deg) # rad, dim (sat, )
    AOP =  np.deg2rad(AOP_deg) # rad, dim (sat, )
    MA =  np.deg2rad(MA_deg) # rad, dim (sat, )

    E0_rad, TA0_rad = E_and_TA_from_MA(MA, e) # rad, dim (sat, )

    w_rot_T = 2*np.pi/(24*3600) # rad/s

    E = np.zeros( ( len(time_span), len(MA)) ) # dim (t, sat )
    for k, (ecc,sma) in enumerate(zip(e,a)):
        for index, t in enumerate(time_span):
            # Hallo la anomalía excéntrica (rad) para cada instante de tiempo y para cada satélite
            def f(x):
                return x - ecc*np.sin(x) - np.sqrt(mu/sma**3)*t
            E[index, k] = fsolve(f, 0) # rad
    E_rad = np.where(E<0, E+2*np.pi, E)
    E_rad += E0_rad

    # Hallo el radio (km) y la anomalía verdadera (rad)
    r = a*(1-e*np.cos(E_rad)) # km, dim (t, sat)
    TA_rad = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2)) # rad, dim (t, sat)
    TA_rad = np.where(TA_rad<0, TA_rad+2*np.pi, TA_rad)
    TA_rad += TA0_rad

    # Hallo el tiempo sidereo local de la estación 
    # Calculo el LST de la GS en el instante inicial
    time = time_in_constellation_epoch(constellation) # dim (s, )
    LST_deg = np.zeros(len(constellation))
    for index, sat_epoch in enumerate(time):
        LST_deg[index] = sat_epoch.sidereal_time('mean', longitude=long_GS*u.deg).value # h, min, s; dim (s, )
    # Calculo el LST de la GS para todos los instantes temporales
    # Para hacer la suma tengo que hacer broadcasting. Quiero acabar con dimensiones (t, sat)
    # Necesito que la parte del LST_deg sea dim (1, sat) y la parte de time_span sea dim (t, 1)
    # (1, sat) + (t, 1) = (t, sat) 
    LST_rad = LST_deg.reshape(1,len(LST_deg))*15*np.pi/180 + w_rot_T*time_span.reshape(len(time_span),1) 

    # Ec(24) apartado 6.4 Elices, dim r_dot_cospsi (t, sat)
    lat_GS = np.deg2rad(lat_GS)
    r_dot_cospsi = r*( np.cos(lat_GS)*( np.cos(AOP+TA_rad)*np.cos(LST_rad-RAAN) + np.cos(i)*np.sin(AOP+TA_rad)*np.sin(LST_rad-RAAN) ) + np.sin(lat_GS)*np.sin(i)*np.sin(AOP+TA_rad) )

    # Ec(21) apartado 6.4 Elices 
    eps_sat = np.arcsin((r_dot_cospsi - Rt)/np.sqrt(r**2 + Rt**2 - 2*Rt*r_dot_cospsi)) # dim (t, sat)
    eps_sat = np.where(eps_sat<0, eps_sat + 2*np.pi, eps_sat)

    return eps_sat

def GMAT_parameters(constellation):
    RAAN_deg, i_deg, e, a, AOP_deg, MA_deg = from_TLE_to_OrbParams(constellation)    # dim (sat, )

    _, TA0_rad = E_and_TA_from_MA(np.deg2rad(MA_deg), e) # rad, dim (sat, )

    times = time_in_constellation_epoch(constellation)
    
    data = []
    for j in range(len(constellation)):
        data.append([constellation[j], a[j], e[j], i_deg[j], RAAN_deg[j], AOP_deg[j], np.rad2deg(TA0_rad[j]), times[j] ])
    
    tabla = pd.DataFrame(data, columns=["Satelite", "SMA [km]", "ECC", "INC [deg]", "RAAN [deg]", "AOP [deg]", "TA0 [deg]", "Epoch in UTC Gregorian"])
    print(tabla)
    tabla.to_csv('outputs/GMAT_parameters.csv')

def contact_locator(eps_sat, eps_GS, time_span,constellation):

    eps_GS= np.deg2rad(eps_GS) # rad
    status = np.where((eps_sat >= eps_GS) & (eps_sat <= np.pi- eps_GS), 1, 0) # dim (t, sat)

    contact = np.where(status == 1)[0]
    contact_locator = []
    contact_locator.append(contact[0])
    for j in range(2,len(contact)):
        if contact[j]-contact[j-1] > 1:
            contact_locator.append(contact[j-1])
            contact_locator.append(contact[j])
    contact_locator.append(contact[-1])

    start_time = time_in_sat_epoch(constellation)
    start_time.format = 'isot'
    
    duration_vec = []
    with open(f'outputs/{constellation}.csv', 'w') as csvfile:
        fieldnames = ['Start time', 'Stop time', 'Duration [s]']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
        for j in range(0,len(contact_locator),2):
            d = {
                "Start time": start_time + time_span[contact_locator[j]]*u.s, 
                "Stop time": start_time + time_span[contact_locator[j+1]]*u.s,
                "Duration [s]": time_span[contact_locator[j+1]]-time_span[contact_locator[j]]
             }
            writer.writerow(d)
            duration_vec.append(time_span[contact_locator[j+1]]-time_span[contact_locator[j]])

    print("Para el satélite", constellation, ": start time =", start_time)            
    print("Number of contacts =", len(contact_locator)/2)
    print("Duración de los contactos [s]:", duration_vec)
    print("-----------------------------------------------------")