from matplotlib import pyplot as plt
import numpy as np

file_in = "datos_energia.dat"
file_out = "grafica_energia"

datos=np.loadtxt(file_in)

T=datos[:,0]
V=datos[:,1]
E=datos[:,2]

n=len(datos)
x=range(n)

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(x,T,label="Energía cinética")
ax.plot(x,V,label="Energía potencial")
ax.plot(x,E,label="Energía total")
ax.legend()
fig.savefig(file_out)