from generating_FOSSA_TLE_function import generate_FOSSA_TLE
from FoM_functions import elevacion_sat
# from read_celestrak import _generate_url_no_xml
from urllib.request import urlopen
import numpy as np
from generating_TLE_function import *

# La GS está en Reyjkiavik (65.64737 N -20.24609 E) y se considera una elevación mínima de 15 deg
eps_GS = 15 # deg
lat_GS = 65.64737 # deg
long_GS = -20.24609 # deg

# Creo un vector de tiempos (en s) que va de 0 (paso por el perigeo) hasta t_fin en incrementos de dt
t_fin = 23*3600 # s
dt = 2 # s
time_span = np.linspace(0, t_fin, int(t_fin/dt)) # dim (t, )

# FOSSA CONSTELATION --> REAL TLEs
# The satellites catalog numbers are:
catalog_numbers = (
    50984,
    50985,
    50987,
    50990,
    50991,
    51001,
    52750,
    52773,
    52776,
    52778,
    52779
)

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

# real_TLE_list = []
# for index, num in enumerate(catalog_numbers):
#     url = _generate_url_no_xml(catalog_number=num, international_designator=None,name=None)

#     # Abrir URL.
#     r = urlopen(url)
#     # Leer el contenido y e imprimir su tamaño.
#     url_TLE = r.read()
#     # Cerrar para liberar recursos.
#     r.close()

#     encoding = 'utf-8'
#     s = str(url_TLE, encoding)

#     _, TLE1, TLE2 = s.splitlines()

#     TLE_sat = [TLE1, TLE2]

#     real_TLE_list.append(TLE_sat)

# print("real TLE list")
# print(real_TLE_list)

# NOTIONAL SATELLITES --> NOTIONAL TLEs
a_notional = np.array([6852.353270, 6851.270166])
e_notional = np.array([0.001066, 0.001104])
i_notional = np.array([97.4682, 97.4687])
RAAN_notional = np.array([88.7588, 88.6583])
AOP_notional = np.array([26.0355, 28.3519])
E_notional = np.array([334.088644, 331.772313])
MA_notional = E_notional + e_notional*np.sin(E_notional)
ha_notional = a_notional*(1+e_notional)

# _, _, max_t_index, _ = elevacion_sat(time_span, lat_GS, long_GS, constellation)

# accedo al epoch time del satélite que se ha actualizado más recientmente
# epoch_time = real_TLE_list[max_t_index][0].split()[3]
epoch_time = 23016.735540665

Satellite_catalog_number = 52780

notional_TLE_list = []
for sat in range (len(a_notional)):
    notional_TLE_list.append( generate_FOSSA_TLE(i_notional[sat], RAAN_notional[sat], e_notional[sat], AOP_notional[sat], MA_notional[sat], a_notional[sat], epoch_time, str(Satellite_catalog_number)) )
    Satellite_catalog_number += 1

print("notional_TLE_list")
print(notional_TLE_list)