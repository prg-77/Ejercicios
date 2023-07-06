from matplotlib import pyplot as plt
import numpy as np

file_in = "grafica_pendulo.dat"
file_out = "0,001phivspsi-E=15"

datos=np.loadtxt(file_in)

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(datos[:,0], datos[:,2])
ax.set_xlabel(r"$\phi$ (rad)")
ax.set_ylabel(r"$\psi$ (rad)")
#ax.set_xlim(-np.pi/2,np.pi/2)
#ax.set_ylim(-np.pi/2,np.pi/2)

fig.savefig(file_out)