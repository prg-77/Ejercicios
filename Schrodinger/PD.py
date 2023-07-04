from matplotlib import pyplot as plt
import numpy as np

file_in = "PD.dat"
file_out = "grafica_PD"

datos=np.loadtxt(file_in)

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(datos[:,0],datos[:,1])

fig.savefig(file_out)
