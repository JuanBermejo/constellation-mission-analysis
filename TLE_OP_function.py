import numpy as np

from poliastro.core.angles import M_to_E, E_to_nu
from astropy.time import Time
from read_celestrak import load_gp_from_celestrak, sat_to_dict
from scipy.optimize import fsolve

mu = 3.986E5

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
    # mu = 3.986E5
    
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
    """
    Transforms Mean Anomaly into Excentric Anomaly and True Anomaly.

    Parameters
    ----------
    MA_rad : float array
        Mean Anomaly [rad]
    e : float array
        Eccentricity [-]

    Returns
    -------
    E_rad : float array
        Excentric Anomaly [rad]
    TA_rad: float array
        True Anomaly [rad]
    """

    E_rad = np.zeros_like(MA_rad)
    TA_rad = np.zeros_like(MA_rad)

    for index in range(len(MA_rad)):
        E_rad[index] = M_to_E(MA_rad[index], e[index]) # dim (sat, )
        TA_rad[index] = E_to_nu(E_rad[index], e[index]) # dim (sat, )
    TA_rad = np.where(TA_rad <0 , TA_rad + 2*np.pi, TA_rad)

    return E_rad, TA_rad

def time_in_sat_epoch(sat_name):
    """
    Gives the epoch time in UTC of the last update of a determined satellite

    Parameters
    ----------
    sat_name : string
        Name of the sat for which the time is desired

    Returns
    -------
    time : time object
        Epoch time for the satellite in UTC
    """
    sat = list(load_gp_from_celestrak(name=sat_name))[0]
    sat_dict = sat_to_dict(sat,sat_name)
    time = Time(sat_dict["EPOCH"], scale='utc')

    return time

def time_in_constellation_epoch(constellation):
    """
    Gives the epoch time in UTC of the last update of each satellite in a constellation

    Parameters
    ----------
    constellation : string array
        Name of the satellites in the constellation for which the times are desired

    Returns
    -------
    time : list of time objects
        Epoch time for the satellite in UTC
    """

    time = []
    for sat_name in constellation:
        sat = list(load_gp_from_celestrak(name=sat_name))[0]
        sat_dict = sat_to_dict(sat,sat_name)
        time.append(Time(sat_dict["EPOCH"], scale='utc'))

    return time

def prop_to_start_time(constellation, t0, T):
    """
    Given a constellation, calculates the most recent update of the TLE in the constellation and gives the parameters nedeed to propagate
    the rest of the satellites to that time.

    Parameters
    ----------
    constellation : string array
        Name of the satellites in the constellation for which the times are desired
    t0: float array
        Time elapsed [s] from perigee passage to TLE update time for each satellite in a constellation
    T: float array
        Orbital period for each satellite in a constellation [s]
    Returns
    -------
    t_in_orbit : float array
        Time elapsed [s] from perigee passage to start contact evaluation time for each satellite in a constellation
    max_t_index : integer
        Index of the satellite in the constellation which has the most recent update of the TLE
    span_difference_array : float array   
        Propagation time needed for each satellite to reach the epoch time of the satellite which has the most recent update of the TLE
    """

    constellation_times = time_in_constellation_epoch(constellation)
    span_difference = []
    for time in constellation_times:
        difference = max(constellation_times) - time
        span_difference.append(difference.value) # time in days 

    max_t_index = constellation_times.index(max(constellation_times))

    span_difference_array = np.array(span_difference)*24*3600
    t_in_orbit = (t0 + span_difference_array)%T

    return t_in_orbit, max_t_index, span_difference_array

def from_t_to_E_TA(e, a, t_in_orbit):
    """
    For a given time and a determined orbit, gives the excentric anomaly and the true anomaly at that time

    Parameters
    ----------
    a: float array
        Semi-major axis [km]
    e: float array
        Eccenticity [-]
    t_in_orbit : float array
        Time elapsed [s] from perigee passage to start contact evaluation time for each satellite in a constellation

    Returns
    -------
    E_in_orbit : float
        Excentric Anomaly [rad]
    TA_in_orbit : float 
        True Anomaly [rad]
        
    """

    E_in_orbit = np.zeros_like(e) # dim (sat, )
    for k, (ecc,sma) in enumerate(zip(e,a)):
        t = t_in_orbit[k]
        # Hallo la anomalía excéntrica (rad) para cada instante de tiempo y para cada satélite
        def f(x):
            return x - ecc*np.sin(x) - np.sqrt(mu/sma**3)*t
        E_in_orbit[k] = fsolve(f, 0) # rad

    TA_in_orbit = np.zeros_like(E_in_orbit)
    for index in range(len(E_in_orbit)):
        TA_in_orbit[index] = E_to_nu(E_in_orbit[index], e[index]) # dim (sat, )
    TA_in_orbit = np.where(TA_in_orbit <0 , TA_in_orbit + 2*np.pi, TA_in_orbit)

    return E_in_orbit, TA_in_orbit