import numpy as np
import matplotlib.pyplot as plt

def print_FoM(datos_graf_dur, datos_graf_rev, labels_vec):

    # FIGURA DE LA DURACIÓN DE LOS CONTACTOS
    fig_dur = plt.figure(1, figsize=(9, 5))

    # Creación el subgrafico
    ax_dur = fig_dur.add_subplot(111)

    # título del gráfico y label del eje y
    plt.ylabel('Tiempo de visita [min]')
    ax_dur.set_title("Tiempo de visita de los contactos ordenados", loc='center', fontdict={'fontsize':14, 'fontweight':'bold'})
    
    # crear el grafico de cajas, incluir labels de los datos que se están representando
    bp_dur = ax_dur.boxplot(datos_graf_dur, labels = labels_vec)

    # visualización mas facil los puntos atípicos
    for flier in bp_dur['fliers']:
        flier.set(marker='o', color='red', alpha=0.5)

    # GRAFICO DEL TIEMPO DE REVISITA
    fig_rev = plt.figure(2, figsize=(9, 6))

    # Creación el subgrafico
    ax_rev = fig_rev.add_subplot(111)

    # título del gráfico y label del eje y
    plt.ylabel('Tiempo de revisita [min]')
    ax_rev.set_title("Tiempo de revisita de los contactos ordenados", loc='center', fontdict={'fontsize':14, 'fontweight':'bold'})
    
    # crear el grafico de cajas, incluir labels de los datos que se están representando
    bp_rev = ax_rev.boxplot(datos_graf_rev, labels = labels_vec)
    
    # visualización mas facil los puntos atípicos
    for flier in bp_rev['fliers']:
        flier.set(marker='o', color='red', alpha=0.5)

    return plt.show() 
    
    
    
# n_fossa_visits =  np.arange(len(real_contact_table))
# n_fossa_revisits = np.arange(len(real_contact_table)-1)
# n_result_visits =  np.arange(len(result_contact_table))
# n_result_revisits = np.arange(len(result_contact_table)-1)
# fig, ax = plt.subplots()
# ax.plot(n_fossa_visits, real_sort_contacts["Duration [s]"]/60, marker = 'o', color = 'tab:purple')
# ax.plot(n_result_visits, result_sort_contacts["Duration [s]"]/60, marker = 'o', color = 'tab:green', linestyle = 'dotted')
# plt.xlabel('Nº de contacto')
# plt.ylabel('Tiempo de visita [min]')
# ax.set_title("Tiempo de visita de los contactos ordenados", loc='center', fontdict={'fontsize':14, 'fontweight':'bold'})
# # plt.show()
# fig, ax = plt.subplots()
# ax.plot(n_fossa_revisits, real_t_revisita, marker = 'o', color = 'tab:purple')
# ax.plot(n_result_revisits, result_t_revisita, marker = 'o', color = 'tab:green', linestyle = 'dotted')
# plt.xlabel('Nº de contacto')
# plt.ylabel('Tiempo de revisita [min]')
# ax.set_title("Tiempo de revisita entre contactos", loc='center', fontdict={'fontsize':14, 'fontweight':'bold'})
# plt.show()