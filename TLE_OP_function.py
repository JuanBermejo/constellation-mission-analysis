import numpy as np

from astropy.coordinates import CartesianRepresentation, CartesianDifferential
from poliastro.util import time_range
from astropy import units as u
from astropy.time import Time

from read_celestrak import *


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

    return RAAN, i, e, a, AOP