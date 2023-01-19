import pandas as pd
import numpy as np

from astropy.time import Time
from generating_FOSSA_TLE_function import generate_FOSSA_TLE
from FoM_functions import elevacion_sat
from read_celestrak import _generate_url_no_xml, sat_to_dict, load_gp_from_celestrak
from urllib.request import urlopen
from generating_TLE_function import *
from main_pases import get_az_sequence, constellation_TLE_list
from FoM_print import print_FoM

print("start running")
# La GS está en Reyjkiavik (65.64737 N -20.24609 E) y se considera una elevación mínima de 15 deg
eps_GS = 15 # deg
lat_GS = 65.64737 # deg
long_GS = -20.24609 # deg
h = 105  # m
min_el = 15  # deg
GS = [long_GS, lat_GS, h, min_el]

# Creo un vector de tiempos (en s) que va de 0 (paso por el perigeo) hasta t_fin en incrementos de dt
t_fin = 23*3600 # s
dt = 2 # s
time_span = np.linspace(0, t_fin, int(t_fin/dt)) # dim (t, )

# FOSSA CONSTELATION --> REAL TLEs
# The satellites catalog numbers are:
catalog_numbers = ( 50984, 50985, 50987, 50990, 50991, 51001, 52750, 52773, 52776, 52778, 52779 )

constellation = (
    "FOSSASAT-2E2", "FOSSASAT-2E1", "FOSSASAT-2E4", "FOSSASAT-2E3", "FOSSASAT-2E6", "FOSSASAT-2E5",
    "FOSSASAT-2E11", "FOSSASAT-2E12", "FOSSASAT-2E13", "FOSSASAT-2E7", "FOSSASAT-2E8"
)

## CONTACTOS DE SATÉLITES FOSSA
real_TLE_list = constellation_TLE_list(catalog_numbers)
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


# NOTIONAL SATELLITES --> NOTIONAL TLEs
a_notional = np.array([6902.741031727869, 6902.741031727869, 6902.741031727869, 6902.741031727869])
e_notional = np.array([0.0010708090909090908, 0.0010708090909090908, 0.0010708090909090908, 0.0010708090909090908])
i_notional = np.array([97.49506363636364, 97.49506363636364, 97.49506363636364, 97.49506363636364])
RAAN_notional = np.array([0, 0, 0, 0])
AOP_notional = np.array([0, 0, 0, 0])
E_notional = np.array([0, 90, 180, 270])
MA_notional = E_notional + e_notional*np.sin(E_notional)
ha_notional = a_notional*(1+e_notional)
notional_names = ("52780", "52781", "52782", "52783")

_, _, max_t_index, _ = elevacion_sat(time_span, lat_GS, long_GS, constellation)

# accedo al epoch time del satélite que se ha actualizado más recientmente
epoch_time = float(real_TLE_list[max_t_index][0].split()[3])

Satellite_catalog_number = 52780

notional_TLE_list = []
for sat in range (len(a_notional)):
    notional_TLE_list.append( generate_FOSSA_TLE(i_notional[sat], RAAN_notional[sat], e_notional[sat], AOP_notional[sat], MA_notional[sat], a_notional[sat], epoch_time, str(Satellite_catalog_number)) )
    Satellite_catalog_number += 1

## TABLAS DE CONTACTOS
# 1) Hallar el inicio del estudio de los contactos
epoch_sat =  list(load_gp_from_celestrak(name=constellation[max_t_index]))[0]
date = sat_to_dict(epoch_sat, constellation[max_t_index])['EPOCH']
date = Time(date, format='isot', scale='utc')
date.format = 'iso'

# 2) satélites reales de la constelación
df_list = []
for index, TLE in enumerate(real_TLE_list):
    df = get_az_sequence(TLE, GS, str(date), 1, constellation[index])
    df_list.append(df)

real_contact_table = pd.concat(df_list)

# 3) satélites ficticios 
notional_list = []
for index, TLE in enumerate(notional_TLE_list):
    df = get_az_sequence(TLE, GS, str(date), 1, notional_names[index])
    notional_list.append(df)

result_contact_table = pd.concat([real_contact_table, pd.concat(notional_list)])

## DURACIÓN Y TIEMPO DE REVISITA
# Fossa
real_sort_contacts = real_contact_table.sort_values('Start Time [s]')
real_t_revisita = np.diff(real_sort_contacts['Start Time [s]'])
real_t_revisita = real_t_revisita*24*60
real_t_revisita = np.array([t.value for t in real_t_revisita])

# Notional
result_sort_contacts = result_contact_table.sort_values('Start Time [s]')
result_t_revisita = np.diff(result_sort_contacts['Start Time [s]'])
result_t_revisita = result_t_revisita*24*60
result_t_revisita = np.array([t.value for t in result_t_revisita])

## GRÁFICO DE CAJAS
datos_dur_fossa = real_sort_contacts["Duration [s]"]/60
datos_dur_comb = result_sort_contacts["Duration [s]"]/60

datos_rev_fossa = real_t_revisita
datos_rev_comb = result_t_revisita

datos_graf_dur = [datos_dur_fossa, datos_dur_comb]
datos_graf_rev = [datos_rev_fossa, datos_rev_comb]

print_FoM(datos_graf_dur, datos_graf_rev)