import numpy as np
import matplotlib.pyplot as plt

# This files contain the TLEs for each FOSSA satellite between de 1st and the 31th of October 2022
FOSSA_TLE_list = (
    "sat000050984.txt", 
    "sat000050985.txt",
    "sat000050987.txt", 
    "sat000050990.txt",  
    "sat000051001.txt", 
    "sat000052750.txt",	
    "sat000052773.txt",	
    "sat000052776.txt",	
    "sat000052778.txt", 
    "sat000052779.txt",
    # "sat000050991.txt" 
)

# The desired information (First and second time derivative of mean motion + drag term) are gathered in the first TLE line
TLE_1 = np.arange(1, 149, 3)

MM_Dot = np.zeros((len(FOSSA_TLE_list), len(TLE_1)))
MM_DDot = np.zeros((len(FOSSA_TLE_list), len(TLE_1)))
BSTAR = np.zeros((len(FOSSA_TLE_list), len(TLE_1)))
sat = []

for j, text_file in enumerate(FOSSA_TLE_list):
    name = "Historico 1-31 oct FOSSA/" + text_file
    file = open(name, "r") 
    
    content = file.readlines() 
    sat.append(content[0])

    counter = 0
    for i in TLE_1:
        line =  content[i]
        MM_Dot[j, counter] = line[33:42]
        MM_DDOt_string = line[44:50] + "e" + line[50:52]
        MM_DDot[j, counter] = float(MM_DDOt_string)
        BSTAR_string = "0." + line[54:59] +"e" + line[59:61]
        BSTAR[j, counter] = float(BSTAR_string)
        counter += 1

    file.close()

lecturas = np.arange(0, len(TLE_1))
sat_colors = ("tab:blue", "tab:green", "tab:orange", "tab:purple", "tab:red", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan")

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
for i in range(len(FOSSA_TLE_list)):
    ax1.plot(lecturas, MM_Dot[i], color = sat_colors[i], label = sat[i])
    ax2.plot(lecturas, BSTAR[i], color = sat_colors[i], label = sat[i] , linestyle = "--")
    
ax1.set_xlabel("Lecture")
ax1.set_ylabel("First time derivative of mean motion [—]")
ax2.set_ylabel("Drag term [--]")
plt.legend()
plt.title("TLE-terms evolution")
plt.show()