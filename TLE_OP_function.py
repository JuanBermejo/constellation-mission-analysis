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
    w: float array
        Argument of the perigee [deg]
    """

    # Definition of the earth gravitational parameter [km^3/s^2]
    mu = 3.986E5
    
    # Definition of the result matrix for the orbital parameters 
    RAAN = np.zeros(len(constellation))
    i = np.zeros(len(constellation))
    e = np.zeros(len(constellation))
    a = np.zeros(len(constellation))
    w = np.zeros(len(constellation))
    
    # The first loop considers all the satellites in the constellation
    for j in range(len(constellation)): 
        sat = list(load_gp_from_celestrak(name=constellation[j]))[0]

        now = Time.now()
        now.jd1, now.jd2

        error, r_vec, v_vec = sat.sgp4(now.jd1, now.jd2)
        assert error == 0
        
        # The angular momemtum vector is calculated with the cross product between r_vecs and v_vecs
        h_vec = np.cross(r_vec,v_vec)
        
        RAAN[j] = np.arcsin(h_vec[0]/np.sqrt(h_vec[0]**2+h_vec[1]**2))*180/np.pi
        h = np.linalg.norm(h_vec)
        i[j] = np.arccos(h_vec[2]/h)*180/np.pi
        
        e_vec = np.cross(v_vec,h_vec)/mu - r_vec/np.linalg.norm(r_vec)
        e[j] = np.linalg.norm(e_vec)
        
        p = h**2/mu
        a[j] = p/(1-e[j]**2)
        n = np.array( [np.cos(RAAN[j]), np.sin(RAAN[j]), 0] )
        w[j] = np.arccos( np.dot(e_vec,n)/e[j] )*180/np.pi   

    return RAAN, i, e, a, w
        
        