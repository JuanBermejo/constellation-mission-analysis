import numpy as np
import pandas as pd
from astropy import units as u
import csv
from poliastro.core.angles import E_to_nu

from TLE_OP_function import from_TLE_to_OrbParams, E_and_TA_from_MA, time_in_sat_epoch, time_in_constellation_epoch 
from TLE_OP_function import prop_to_start_time, from_t_to_E_TA
from scipy.optimize import fsolve

# DATOS
Rt = 6378 # km
mu = 3.986e5 # km3/s2

def elevacion_sat(time_span, lat_GS, long_GS, constellation):
    """
    Cálculo de la elevación en cada satélite en una constelación para unos instantes temporales y una GS

    Parámeteros
    ----------
    time_span : float array
        instantes temporales de evaluación de los contactos
    lat_GS: float
        latitud de la Ground Station [deg]
    constellation : list of strings
        names of the satellites in the constellation

    Devuelve
    -------
    eps_sat : float array
        Elevación del satélite en cada instante temporal [rad]
    """

    # Obtengo los parámetros orbitales para la última actualización del TLE de cada satélite
    RAAN_deg, i_deg, e, a, AOP_deg, MA_deg = from_TLE_to_OrbParams(constellation)    # dim (sat, )

    # Paso los datos de los parámetros orbitales a radianes para utilizar np.cos, np.sin, np.tan, etc.
    RAAN =  np.deg2rad(RAAN_deg) # rad, dim (sat, )
    i =  np.deg2rad(i_deg) # rad, dim (sat, )
    AOP =  np.deg2rad(AOP_deg) # rad, dim (sat, )
    MA =  np.deg2rad(MA_deg) # rad, dim (sat, )

    # Para el instante de la actualización de cada TLE, calculo la anomalía excéntrica (E0) y el tiempo respecto del paso 
    # por el perigeo (t0) con la ecuación de Kepler. Cálculo del periodo orbital (T) de cada satélite
    E0_rad, TA0_rad = E_and_TA_from_MA(MA, e) # rad, dim (sat, )
    t0 = (E0_rad-e*np.sin(E0_rad))*np.sqrt(a**3/mu)
    T = 2*np.pi*np.sqrt(a**3/mu)

    # Hasta aquí tenemos unos satélites de los que conocemos sus parámetros orbitales a través de sus TLEs.
    # Los TLEs no se actualizan al mismo tiempo, así que hallamos el satélite que ha sido actualizado hace menos tiempo y vemos 
    # cuánto tiempo falta por propagar al resto para llegar hasta ese instante temporal en segundos (delta_t).
    # Con el dato de la MA del TLE sacamos la E del instante del TLE (E0) y hallamos el tiempo para el que sucede (t0)
    # Con el instante inicial y lo que falta para llegar al momento temporal del último satélite actualizado hallamos donde está
    # el satélite en el momento en el que vamos a empezar a evaluar los contactos (t_in_orbit).
    # Para ese instante temporal obtenemos la E (E_in_orbit).
    # Se puede comparar con plot_TLE_OP.ipynb si la distribución de satélites es coherente.  
      
    t_in_orbit, max_t_index, delta_t = prop_to_start_time(constellation, t0, T)
    E_in_orbit, TA_in_orbit = from_t_to_E_TA(e, a, t_in_orbit)

    E_in_orbit = np.zeros_like(e) # dim (sat, )

    w_rot_T = 2*np.pi/(24*3600) # rad/s

    # Ahora se va a realizar la propagación de los satélites con la ecuación de Kepler para evaluar los contactos. 
    # Para ello, para cada satélite se halla la anomalía excéntrica en unos instantes temporales (time_span) desde 0 hasta el
    # tiempo que se halla definido en FoM_main.py . Como el tiempo t = 0 s representa el paso del satélite por el perigeo, se
    # debe sumar la anomalía excéntrica a la que se encontraba el satélite antes de empezar la evaluación de los contactos 
    # (E_in_orbit). El resultado de la propagación con la ecuación de Kepler (E) se emplea, por un lado, para saber la anomalía
    # excéntrica en cada punto de la propagación (E_rad = E + E_in_orbit) y, por otro, para hallar la anomalía verdadera (TA_rad)

    E = np.zeros( ( len(time_span), len(MA)) ) # dim (t, sat )
    for k, (ecc,sma) in enumerate(zip(e,a)):
        for index, t in enumerate(time_span):
            # Hallo la anomalía excéntrica (rad) para cada instante de tiempo y para cada satélite
            def f(x):
                return x - ecc*np.sin(x) - np.sqrt(mu/sma**3)*t
            E[index, k] = fsolve(f, 0) # rad
    E_rad = np.where(E<0, E+2*np.pi, E)
    E_rad += E_in_orbit

    # Hallo el radio (km) y la anomalía verdadera (rad) en cada instante temporal para cada satélite
    r = a*(1-e*np.cos(E_rad)) # km, dim (t, sat)
    TA_rad = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2)) # rad, dim (t, sat)
    TA_rad = np.where(TA_rad<0, TA_rad+2*np.pi, TA_rad)

    TA_rad += TA_in_orbit 

    # Hallo el tiempo sidereo local de la estación en el inicio de la evaluación de los contactos --> en el momento 
    # temporal del satélite que ha actualizado hace menos tiempo su TLE (GS_time). GS_time es un objeto de astropy.time que tiene 
    # asociado un atributo que convierte el tiempo dado a TLS para una longitud determinada ( .sideral_time() ). El resultado de 
    # este atributo es un tiempo en horas, minutos y segundos que para el problema que se está evaluando hay que pasar a radianes.
    # Además, se necesita un LST de la GS para cada instante de evaluación de los contactos --> GS_LST_rad
    GS_time = time_in_sat_epoch(constellation[max_t_index])
    GS_LST_deg = GS_time.sidereal_time('mean', longitude=long_GS*u.deg).value # h, min, s; 
    GS_LST_rad = GS_LST_deg*15*np.pi/180 + w_rot_T*time_span # dim (t, )

    # Con todos estos datos se calcula la elevación del satélite 
    # Ec(24) apartado 6.4 Elices. dim r_dot_cospsi (t, sat)
    lat_GS = np.deg2rad(lat_GS)
    GS_LST_rad = GS_LST_rad.reshape(len(GS_LST_rad),1) # dim (t, 1)
    r_dot_cospsi = r*( np.cos(lat_GS)*( np.cos(AOP+TA_rad)*np.cos(GS_LST_rad-RAAN) + np.cos(i)*np.sin(AOP+TA_rad)*np.sin(GS_LST_rad-RAAN) ) + np.sin(lat_GS)*np.sin(i)*np.sin(AOP+TA_rad) )
    # Ec(21) apartado 6.4 Elices 
    eps_sat = np.arcsin((r_dot_cospsi - Rt)/np.sqrt(r**2 + Rt**2 - 2*Rt*r_dot_cospsi)) # dim (t, sat)
    eps_sat = np.where(eps_sat<0, eps_sat + 2*np.pi, eps_sat)

    return eps_sat, delta_t

def GMAT_parameters(constellation):
    """
    Print de una tabla con los parámetros orbitales de cada satélite en la constelación y su epoch en UTC Gregorian

    Parámeteros
    ----------
    constellation : list of strings
        names of the satellites in the constellation
    """
    
    RAAN_deg, i_deg, e, a, AOP_deg, MA_deg = from_TLE_to_OrbParams(constellation)    # dim (sat, )

    _, TA0_rad = E_and_TA_from_MA(np.deg2rad(MA_deg), e) # rad, dim (sat, )

    TLE_times = time_in_constellation_epoch(constellation)
    
    data = []
    for j in range(len(constellation)):
        data.append([constellation[j], a[j], e[j], i_deg[j], RAAN_deg[j], AOP_deg[j], np.rad2deg(TA0_rad[j]), TLE_times[j] ])
    
    tabla = pd.DataFrame(data, columns=["Satelite", "SMA [km]", "ECC", "INC [deg]", "RAAN [deg]", "AOP [deg]", "TA0 [deg]", "Epoch in UTC Gregorian"])
    print(tabla)
    # tabla.to_csv('outputs/GMAT_parameters.csv')

def contact_locator(eps_sat, eps_GS, time_span, sat_name, delta_t):
    """
    Print de una tabla con la información de los contactos para un satélite dado

    Parámeteros
    ----------
    eps_sat : float array
        elevación del satélite en cada instante temporal de evaluación de los contactos
    eps_GS: float
        elevación mínima de la Ground Station [deg]
    time_span : float array
        instantes temporales de evaluación de los contactos
    sat_name: string
        nombre del satélite del que se están evaluando los contactos
    delta_t: float
        desfase temporal desde la actualización del TLE del satélite hasta el inicio de la evaluación de los contactos
    """
    
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

    start_time = time_in_sat_epoch(sat_name) + delta_t*u.s
    start_time.format = 'isot'
    
    data = []
    for j in range(0,len(contact_locator),2):
        data.append([ sat_name, start_time + time_span[contact_locator[j]]*u.s, start_time + time_span[contact_locator[j+1]]*u.s, time_span[contact_locator[j+1]]-time_span[contact_locator[j]] ])
    tabla = pd.DataFrame(data, columns=["Satellite", "Start Time [s]", "Stop Time [s]", "Duration [s]"])
    print(tabla)
    print("-----------------------------------------------------")