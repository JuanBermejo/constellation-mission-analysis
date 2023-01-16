import numpy as np
from FoM_functions import elevacion_sat, GMAT_parameters, contact_locator, contact_locator_notional_sat
from TLE_OP_function import time_in_sat_epoch
import matplotlib.pyplot as plt
import pandas as pd

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

eps_sat, delta_t, max_t_index, TA_in_orbit = elevacion_sat(time_span, lat_GS, long_GS, constellation)

print("Los parámetros orbitales obtenidos de los TLEs de los satélites son:")
GMAT_parameters(constellation)

## CONTACTOS DE SATÉLITES FOSSA
contact_table = contact_locator(eps_sat, eps_GS, time_span, constellation, delta_t)
# print("Los contactos de cada satélite FOSSA son los siguientes:")
# print(contact_table)
# print("----------------------------------------------------------------------------")
sort_contacts = contact_table.sort_values('Start Time [s]')
print("Ordenando los contactos por cuándo empiezan, se tiene la siguiente secuencia")
print("----------------------------------------------------------------------------")
print(sort_contacts)
t_revisita = np.diff(sort_contacts['Start Time [s]'])
t_revisita = t_revisita*24*60
t_revisita = np.array([t.value for t in t_revisita])
# print("El tiempo de revisita que se obtiene es, por lo tanto, el siguiente:")
# print(t_revisita)
# print("----------------------------------------------------------------------------")

## CONTACTOS DE SATÉLITES FICTICIOS
a_notional = np.array([6851.670117, 6877.625962])
e_notional = np.array([0.001033, 0.000887])
i_notional = np.array([97.4647, 97.5325])
RAAN_notional = np.array([88.0386, 134.0995])
AOP_notional = np.array([35.0131, 142.3066])
# TA0_notional = np.array([np.rad2deg(5.67425834), np.rad2deg(1.22055784)])
TA0_notional = np.array([np.rad2deg(0), np.rad2deg(0)])

orbital_params = np.array( [ a_notional , e_notional, i_notional, RAAN_notional, AOP_notional, TA0_notional ] )
notional_names = ("FOSSANOT-2E14", "FOSSANOT-2E15")
epoch_time = time_in_sat_epoch(constellation[max_t_index])
notional_contact_table = contact_locator_notional_sat(orbital_params, notional_names, epoch_time, time_span, lat_GS, long_GS, eps_GS)
# print("Los contactos de cada satélite ficticio son los siguientes:")
# print(notional_contact_table)
# print("----------------------------------------------------------------------------")

## CONTACTOS DE SATÉLITES FOSSA + FICTICIOS

contact_result = pd.concat([contact_table, notional_contact_table])
# print("Los contactos de los satélites ficticios y reales son los siguientes:")
# print(contact_result)
# print("----------------------------------------------------------------------------")

sort_contact_result = contact_result.sort_values('Start Time [s]')
print("Los contactos de los satélites ficticios y reales ordenados son los siguientes:")
print(sort_contact_result)
print("----------------------------------------------------------------------------")
t_revisita_result = np.diff(sort_contact_result['Start Time [s]'])
t_revisita_result = t_revisita_result*24*60
t_revisita_result = np.array([t.value for t in t_revisita_result])

n_fossa_visits =  np.arange(len(contact_table))
n_fossa_revisits = np.arange(len(contact_table)-1)

n_result_visits =  np.arange(len(contact_result))
n_result_revisits = np.arange(len(contact_result)-1)

fig, ax = plt.subplots()
ax.plot(n_fossa_visits, sort_contacts["Duration [s]"]/60, marker = 'o', color = 'tab:purple')
ax.plot(n_result_visits, sort_contact_result["Duration [s]"]/60, marker = 'o', color = 'tab:green', linestyle = 'dotted')
plt.xlabel('Nº de contacto')
plt.ylabel('Tiempo de visita [min]')
ax.set_title("Tiempo de visita de los contactos ordenados", loc='center', fontdict={'fontsize':14, 'fontweight':'bold'})
# plt.show()

fig, ax = plt.subplots()
ax.plot(n_fossa_revisits, t_revisita, marker = 'o', color = 'tab:purple')
ax.plot(n_result_revisits, t_revisita_result, marker = 'o', color = 'tab:green', linestyle = 'dotted')
plt.xlabel('Nº de contacto')
plt.ylabel('Tiempo de revisita [min]')
ax.set_title("Tiempo de revisita entre contactos", loc='center', fontdict={'fontsize':14, 'fontweight':'bold'})
plt.show()