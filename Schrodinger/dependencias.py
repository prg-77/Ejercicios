from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
import math as m

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

file_in = "dependencias.dat"
file_out = "grafica_dependencias"

datos=np.genfromtxt(file_in, missing_values='-')

fig=plt.figure(figsize=[30,20])

landa=0.5
n_ciclos=100


ax1=fig.add_subplot(221)
ax1.errorbar(datos[:,0],datos[:,1], yerr=datos[:,2], color="green", ecolor="red", fmt="o-", capsize=2, markersize=1.5, linewidth=0.5, elinewidth=0.3, label="Valor simulado")
x_values = range(500, 2001, 500)

y_values = [4*(1-landa)/(4*(1-landa)+landa*landa*(m.sin(2*m.pi*n_ciclos*m.sqrt(1-landa)/5))**2) for x in x_values]

plt.plot(x_values, y_values, color="blue", linewidth=0.5, label="Valor teórico")
#ax1.plot(datos[:0], 4*(1-landa)/(4*(1-landa)+landa*landa*(m.sin(2*m.pi*n_ciclos*m.sqrt(1-landa)/5))**2))
ax1.set_xlabel("N")
ax1.set_ylabel("Coeficiente de transmisión")
ax1.legend(loc="upper left")

ax2=fig.add_subplot(222)
ax2.errorbar(datos[:,3],datos[:,4], yerr=datos[:,5], color="green", ecolor="red", fmt="o-", capsize=2, markersize=1.5, linewidth=0.5, elinewidth=0.3, label="Valor simulado")
x_values = np.arange(0.1, 1, 0.001)
y_values = [4*(1-x)/(4*(1-x)+x*x*(m.sin(2*m.pi*n_ciclos*m.sqrt(1-x)/5))**2) for x in x_values]
plt.plot(x_values, y_values, color="blue", linewidth=0.5, label="Valor teórico")
ax2.set_xlabel(r'$\lambda$')
ax2.set_ylabel("Coeficiente de transmisión")
ax2.legend(loc="lower left")

ax3=fig.add_subplot(223)
ax3.errorbar(datos[:,6],datos[:,7], yerr=datos[:,8], color="green", ecolor="red", fmt="o-", capsize=2, markersize=1.5, linewidth=0.5, elinewidth=0.3, label="Valor simulado")
x_values = np.arange(500, 2001, 1)

y_values = [4*(1-landa)/(4*(1-landa)+landa*landa*(m.sin(2*m.pi*(x/8)*m.sqrt(1-landa)/5))**2) for x in x_values]

plt.plot(x_values, y_values, color="blue", linewidth=0.5, label="Valor teórico")
#ax1.plot(datos[:0], 4*(1-landa)/(4*(1-landa)+landa*landa*(m.sin(2*m.pi*n_ciclos*m.sqrt(1-landa)/5))**2))
ax3.set_xlabel("N")
ax3.set_ylabel("Coeficiente de transmisión")
ax3.set_ylim(0.6, 1.2)
ax3.legend(loc="upper left")

fig.savefig(file_out)
