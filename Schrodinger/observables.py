from matplotlib import pyplot as plt
import numpy as np

file_in = "graficas_observables.dat"
file_out = "grafica_obsevables"

datos=np.loadtxt(file_in)

fig=plt.figure(figsize=[20,30])

ax1=fig.add_subplot(321)
ax1.errorbar(datos[:,0],datos[:,1], yerr=datos[:,2], color="green", ecolor="red", fmt="o-", capsize=2, markersize=1.5, linewidth=0.5, elinewidth=0.3)
ax1.set_xlabel("Tiempo")
ax1.set_ylabel("Posición")

ax2=fig.add_subplot(322)
ax2.errorbar(datos[:,0],datos[:,3],yerr=datos[:,4], color="green", ecolor="red", fmt="o-", capsize=2, markersize=1.5, linewidth=0.5, elinewidth=0.3)
ax2.set_xlabel("Tiempo")
ax2.set_ylabel("Momento")

ax3=fig.add_subplot(323)
ax3.errorbar(datos[:,0],datos[:,5],yerr=datos[:,6], color="green", ecolor="red", fmt="o-", capsize=2, markersize=1.5, linewidth=0.5, elinewidth=0.3)
ax3.set_xlabel("Tiempo")
ax3.set_ylabel("Energía potencial")

ax4=fig.add_subplot(324)
ax4.errorbar(datos[:,0],datos[:,7], yerr=datos[:,8], color="green", ecolor="red", fmt="o-", capsize=2, markersize=1.5, linewidth=0.5, elinewidth=0.3)
ax4.set_xlabel("Tiempo")
ax4.set_ylabel("Energía cinética")

ax5=fig.add_subplot(325)
ax5.errorbar(datos[:,0],datos[:,9], yerr=datos[:,10], color="green", ecolor="red", fmt="o-", capsize=2, markersize=1.5, linewidth=0.5, elinewidth=0.3)
ax5.set_xlabel("Tiempo")
ax5.set_ylabel("Energía total")


fig.savefig(file_out)