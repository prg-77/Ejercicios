from matplotlib import pyplot as plt
import numpy as np
import math as m

file_in = "nbarreras.dat"
file_out = "nbarreras"

datos=np.loadtxt(file_in)

fig=plt.figure()
ax=fig.add_subplot(111)
ax.errorbar(datos[:,0],datos[:,1], yerr=datos[:,2], color="green", ecolor="red", fmt="o", capsize=2, markersize=1.5, linewidth=0.5, elinewidth=0.3, label="Datos simulados")
ax.set_xlabel("Número de barreras")
ax.set_ylabel("Coeficiente de transmisión")
def func(x):
    return 0.225*m.exp(-0.0331*x)
x_values = datos[:,0]
y_values = [func(x) for x in x_values]
ax.plot(x_values, y_values, color="blue", linewidth=0.5, label="Ajuste exponencial: $0.225e^{-0.0331x}$")
ax.legend()
fig.savefig(file_out)
