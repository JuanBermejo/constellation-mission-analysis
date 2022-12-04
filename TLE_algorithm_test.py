import numpy as np

from astropy.time import Time
from poliastro.util import time_range
from poliastro.bodies import Earth
from poliastro.plotting import OrbitPlotter3D

from read_celestrak import *
from Ephem_from_GP import *
from generating_TLE_function import *

# Definition of the names in the FOSSA FOSSA_list, as presented in celestrak 
FOSSA_list = (
    "FOSSASAT-2E2", 
    "FOSSASAT-2E1",
    "FOSSASAT-2E4", 
    "FOSSASAT-2E3", 
    "FOSSASAT-2E6", 
    "FOSSASAT-2E5", 
    "FOSSASAT-2E11",	
    "FOSSASAT-2E12",	
    "FOSSASAT-2E13",	
    "FOSSASAT-2E7", 
    "FOSSASAT-2E8", 
)

# Definition of the earth gravitational parameter [km^3/s^2]
mu = 3.986004418E5

# Definition of the result matrix for the orbital parameters 
RAAN = np.zeros(len(FOSSA_list))
i = np.zeros(len(FOSSA_list))
e = np.zeros(len(FOSSA_list))
a = np.zeros(len(FOSSA_list))
w = np.zeros(len(FOSSA_list))
MM = np.zeros(len(FOSSA_list))

RAAN_tle = np.zeros(len(FOSSA_list))
i_tle = np.zeros(len(FOSSA_list))
e_tle = np.zeros(len(FOSSA_list))
MM_tle = np.zeros(len(FOSSA_list))
a_tle = np.zeros(len(FOSSA_list))
w_tle = np.zeros(len(FOSSA_list))

for j in range(len(FOSSA_list)): 
    sat = list(load_gp_from_celestrak(name=FOSSA_list[j]))[0]
    sat_dict = sat_to_dict(sat,FOSSA_list[j])
    now = Time(sat_dict["EPOCH"], scale='utc')
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
    MM[j]  = MeanMotion(a[j])

    # n_vec = np.cross([0,0,1], h_vec)
    n_vec = np.array( [np.cos(np.deg2rad(RAAN[j])), np.sin(np.deg2rad(RAAN[j])), 0] )
    n = np.linalg.norm(n_vec)
    w[j] = np.arccos( np.dot(e_vec,n_vec)/(e[j]*n) )*180/np.pi   

    RAAN_tle[j] = sat_dict["RA_OF_ASC_NODE"]
    i_tle[j] = sat_dict["INCLINATION"]
    e_tle[j] = sat_dict["ECCENTRICITY"]
    w_tle[j] = sat_dict["ARG_OF_PERICENTER"]
    MM_tle[j] = sat_dict["MEAN_MOTION"]
    a_tle[j] = (mu/(MM_tle[j]*2*np.pi/86400)**2)**(1/3)

print(w-w_tle)
print(w)
print(w_tle)