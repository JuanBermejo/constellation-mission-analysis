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

# The first and second time derivative of mean motion + drag term are gathered in the first TLE line
TLE_1 = np.arange(1, 149, 3)
# The eccentricity and the mean motion are gathered in the second TLE line
TLE_2 = np.arange(2, 150, 3)

MM_Dot = np.zeros((len(FOSSA_TLE_list), len(TLE_1)))
MM_DDot = np.zeros((len(FOSSA_TLE_list), len(TLE_1)))
BSTAR = np.zeros((len(FOSSA_TLE_list), len(TLE_1)))
MM = np.zeros((len(FOSSA_TLE_list), len(TLE_2)))
ECC = np.zeros((len(FOSSA_TLE_list), len(TLE_2)))
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
        MM_DDot[j, counter] = line[44:50] + "e" + line[50:52]
        BSTAR[j, counter] = "0." + line[54:59] +"e" + line[59:61]
        counter += 1

    counter = 0
    for i in TLE_2:
        line =  content[i]
        ECC[j, counter] = "0." + line[26:34]
        MM[j, counter] = line[52:64]
        counter += 1
        
    file.close()

sat_colors = ("tab:blue", "tab:green", "tab:orange", "tab:purple", "tab:red",
 "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan")

# plot for the evolution of the first derivative mean motion and the drag term
lecturas = np.arange(0, len(TLE_1))
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

# mean values for the first and second mean motion derivatives and the drag term
MM_Dot_mean = np.mean(MM_Dot)
MM_DDot_mean = np.mean(MM_DDot)
BSTAR_mean = np.mean(BSTAR)

print("MM_Dot_mean =", MM_Dot_mean)
print("MM_DDot_mean =", MM_DDot_mean)
print("BSTAR_mean =", BSTAR_mean)

# plot for the eccentricity vs mean motion 
for i in range(len(FOSSA_TLE_list)):
    plt.scatter(MM[i], ECC[i], color = sat_colors[i], label = sat[i])

plt.xlabel("Mean Motion [rev/day]")
plt.ylabel("Eccentricity [deg]")
plt.legend()
plt.title("Mean Motion VS Eccentricity")
plt.show()