import numpy as np
from FoM_functions import elevacion_sat, GMAT_parameters, contact_locator, time_in_sat_epoch

print("Start running")
# la GS está en Reyjkiavik (65.64737 N -20.24609 E) y se considera una elevación mínima de 15 deg
eps_GS = 15 # deg
lat_GS = 65.64737 # deg
long_GS = -20.24609 # deg

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

t_fin = 23*3600 # s
dt = 2 # s
# Creo un vector de tiempos (en s) que va de 0 (paso por el perigeo) hasta t_fin en incrementos de dt
time_span = np.linspace(0, t_fin, int(t_fin/dt)) # dim (t, )

eps_sat, delta_t = elevacion_sat(time_span, lat_GS, long_GS, constellation)
GMAT_parameters(constellation)
contact_table = contact_locator(eps_sat, eps_GS, time_span, constellation, delta_t)
print(contact_table)
