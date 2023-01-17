import numpy as np
from TLE_OP_function import *
from generating_TLE_function import MeanMotion
import matplotlib.pyplot as plt

# Definition of the names in the FOSSA constellation, as presented in celestrak 
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
    "FOSSASAT-2E8" 
)

# The function "from_TLE_to_OrbParams" returns the instant orbital parameters for a given constellation
RAAN, i, e, a, AOP, _ = from_TLE_to_OrbParams(FOSSA_list)

# Check if the orbit is SS
J2 = 1.0827E-3 #[-]
Rt = 6378 #[km]
alfa_dot = 2*np.pi/(24*3600*365.25)
mu = 3.986E5 #[km3/s2]

a_SS = (-3*Rt**2*J2*np.sqrt(mu)*np.cos(np.radians(i))/(2*alfa_dot*(1-e)**2))**(2/7)
print("For FOSSA satellites, a_SS - a =", a_SS-a)

i_SS = np.arccos(-2*alfa_dot*np.sqrt(a**3/mu)*a**2*(1-e**2)**2/(3*J2*Rt**2))*180/np.pi
print("For FOSSA satellites, i_SS - i =", i_SS-i)

plane_rotation = -3*J2*np.sqrt(mu/a**3)*Rt**2*np.cos(np.radians(i))/(2*a**2*(1-e)**2)*180/np.pi*24*3600
print("The rotation of the orbital plane in one day is equal to", plane_rotation)

# Time of pass through the Ecuator 
# The Righ Ascension of the Sun on the 21/11/2022 is 15h 43m 50s = 235.96 deg 
H0 = (RAAN - 235.96)/15 + 12
print("RAAN =", RAAN)
print("The satellite crosses the Ecuator at", 12+H0)