from matplotlib import pyplot as plt
import numpy as np

file_in = "dependencias.dat"
file_out = "grafica_dependencias"

datos=np.genfromtxt(file_in, missing_values='-')

fig=plt.figure(figsize=[20,20])

ax1=fig.add_subplot(221)
ax1.errorbar(datos[:,0],datos[:,1], yerr=datos[:,2], color="green", ecolor="red", fmt="o-", capsize=2, markersize=1.5, linewidth=0.5, elinewidth=0.3)
ax1.set_xlabel("N")
ax1.set_ylabel("Coeficiente de transmisión")

ax2=fig.add_subplot(222)
ax2.errorbar(datos[:,3],datos[:,4], yerr=datos[:,5], color="green", ecolor="red", fmt="o-", capsize=2, markersize=1.5, linewidth=0.5, elinewidth=0.3)
ax2.set_xlabel("lambda")
ax2.set_ylabel("Coeficiente de transmisión")

fig.savefig(file_out)
