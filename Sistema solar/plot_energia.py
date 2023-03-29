from matplotlib import pyplot as plt
import numpy as np

file_in = "datos_energia.txt"
file_out = "grafica_energia"

datos=np.loadtxt(file_in)

n=len(datos)
x=range(n)

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(x,datos)

fig.savefig(file_out)


