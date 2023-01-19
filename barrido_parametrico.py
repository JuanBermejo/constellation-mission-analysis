import numpy as np
from FoM_functions import GMAT_parameters

constellation = (
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

GMAT_params = GMAT_parameters(constellation)
print(GMAT_params)
e = np.mean(GMAT_params["ECC"])
print("e = ", e)
i = np.mean(GMAT_params["INC [deg]"])
print("i =", i)
J2 = 1.0827e-3
mu = 3.986e5
alfa_dot = 2*np.pi/(24*3600*365.25)
Rt = 6378
a = ( -3*J2*np.sqrt(mu)*Rt**2*np.cos( np.deg2rad(i) )/( 2*alfa_dot*(1-e**2)**2 ) )**(2/7)
print("a =", a)