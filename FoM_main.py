import numpy as np
from FoM_functions import elevacion_sat, GMAT_parameters, contact_locator, contact_locator_notional_sat
from TLE_OP_function import time_in_sat_epoch
import pandas as pd
from FoM_print import print_FoM

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

t_fin = 24*3600 # s
dt = 1 # s
# Creo un vector de tiempos (en s) que va de 0 (paso por el perigeo) hasta t_fin en incrementos de dt
time_span = np.linspace(0, t_fin, int(t_fin/dt)) # dim (t, )

eps_sat, delta_t, max_t_index, TA_in_orbit = elevacion_sat(time_span, lat_GS, long_GS, constellation)

# print("Los parámetros orbitales obtenidos de los TLEs de los satélites son:")
# GMAT_params = GMAT_parameters(constellation)
# print("----------------------------------------------------------------------------")

## CONTACTOS DE SATÉLITES FOSSA
real_contact_table = contact_locator(eps_sat, eps_GS, time_span, constellation, delta_t)
real_sort_contacts = real_contact_table.sort_values('Start Time [s]')
real_t_revisita = np.diff(real_sort_contacts['Start Time [s]'])
real_t_revisita = real_t_revisita*24*60
real_t_revisita = np.array([t.value for t in real_t_revisita])

## CONTACTOS DE SATÉLITES FICTICIOS
a_notional = np.array([6902.741031727869, 6902.741031727869, 6902.741031727869, 6902.741031727869])
e_notional = np.array([0.0010708090909090908, 0.0010708090909090908, 0.0010708090909090908, 0.0010708090909090908])
i_notional = np.array([97.49506363636364, 97.49506363636364, 97.49506363636364, 97.49506363636364])
RAAN_notional = np.array([0, 0, 0, 0])
AOP_notional = np.array([0, 0, 0, 0])
TA0_notional = np.array([0, 90, 180, 270])

orbital_params = np.array( [ a_notional , e_notional, i_notional, RAAN_notional, AOP_notional, TA0_notional ] )
notional_names = ("52780", "52781", "52782", "52783")
epoch_time = time_in_sat_epoch(constellation[max_t_index])
notional_contact_table = contact_locator_notional_sat(orbital_params, notional_names, epoch_time, time_span, lat_GS, long_GS, eps_GS)

## CONTACTOS DE SATÉLITES FOSSA + FICTICIOS

result_contact_table = pd.concat([real_contact_table, notional_contact_table])

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

print("Valor medio del tiempo de revisita en la constelación de FOSSA =", np.mean(real_t_revisita), "minutos")
print("Valor medio del tiempo de revisita en la constelación ficticia RAAN 0º =", np.mean(result_t_revisita), "minutos")
print_FoM(datos_graf_dur, datos_graf_rev, ["FOSSA", "RAAN = 0\u00B0"])