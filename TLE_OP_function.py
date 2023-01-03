import numpy as np

from poliastro.core.angles import M_to_E, E_to_nu
from astropy.time import Time
from read_celestrak import load_gp_from_celestrak, sat_to_dict


def from_TLE_to_OrbParams(constellation):
    """
    Calculates the instant orbital parameters for a given list of satellites.

    Parameters
    ----------
    constellation : list of strings
        names of the satellites in the constellation

    Returns
    -------
    RAAN : float array
        Right ascension of the ascending node [deg]
    i: float array
        Inclination [deg]
    e: float array
        Eccenticity [-]
    a: float array
        Semi-major axis [km]
    AOP: float array
        Argument of the perigee [deg]
    MA: float array
        Mean anomaly [deg]
    """

    # Definition of the earth gravitational parameter [km^3/s^2]
    mu = 3.986E5
    
    # Definition of the result matrix for the orbital parameters 
    RAAN = np.zeros(len(constellation))
    i = np.zeros(len(constellation))
    e = np.zeros(len(constellation))
    MM = np.zeros(len(constellation))
    a = np.zeros(len(constellation))
    AOP = np.zeros(len(constellation))
    MA = np.zeros(len(constellation))
    
    # The first loop considers all the satellites in the constellation
    for j in range(len(constellation)): 
        sat = list(load_gp_from_celestrak(name=constellation[j]))[0]
        sat_dict = sat_to_dict(sat,constellation[j]) 

        RAAN[j] = sat_dict["RA_OF_ASC_NODE"]
        i[j] = sat_dict["INCLINATION"]
        e[j] = sat_dict["ECCENTRICITY"]
        AOP[j] = sat_dict["ARG_OF_PERICENTER"]
        MM[j] = sat_dict["MEAN_MOTION"]
        a[j] = (mu/(MM[j]*2*np.pi/86400)**2)**(1/3) 
        MA[j] = sat_dict["MEAN_ANOMALY"]
    # a, e, i_deg, RAAN_deg, AOP_deg, MA_deg = [6857.26216168], [0.0011313], [97.4697], [74.2364], [77.3213] #TODO
    return RAAN, i, e, a, AOP, MA

def E_and_TA_from_MA(MA_rad, e):
    E_rad = np.zeros_like(MA_rad)
    TA_rad = np.zeros_like(MA_rad)
    for index in range(len(MA_rad)):
        E_rad[index] = M_to_E(MA_rad[index], e[index]) # dim (sat, )
        TA_rad[index] = E_to_nu(E_rad[index], e[index]) # dim (sat, )
    TA_rad = np.where(TA_rad <0 , TA_rad + 2*np.pi, TA_rad)

    return E_rad, TA_rad

def time_in_sat_epoch(sat_name):

    sat = list(load_gp_from_celestrak(name=sat_name))[0]
    sat_dict = sat_to_dict(sat,sat_name)
    time = Time(sat_dict["EPOCH"], scale='utc')
    # time =  Time('2023-01-02T15:46:51.437', scale='utc') #TODO
    return time

def time_in_constellation_epoch(constellation):
    time = []
    for sat_name in constellation:
        sat = list(load_gp_from_celestrak(name=sat_name))[0]
        sat_dict = sat_to_dict(sat,sat_name)
        time.append(Time(sat_dict["EPOCH"], scale='utc'))
    # time =  Time('2023-01-02T15:46:51.437', scale='utc') #TODO
    return time