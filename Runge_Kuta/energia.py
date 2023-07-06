from matplotlib import pyplot as plt
import numpy as np

file_in = "energia.dat"
file_out = "grafica_energia"

datos=np.loadtxt(file_in)

n=len(datos)
x=range(n)

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(x,datos)
plt.xlabel("Tiempo")
plt.ylabel("H-w*p_phi")

fig.savefig(file_out)