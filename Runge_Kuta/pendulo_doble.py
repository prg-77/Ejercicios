from matplotlib import pyplot as plt
import numpy as np

file_in = "grafica_pendulo.dat"
file_out = "prueba"

datos=np.loadtxt(file_in)

fig=plt.figure()

#phi vs psi
'''
ax=fig.add_subplot(111)
ax.plot(datos[:,0], datos[:,2])
ax.set_xlabel(r"$\phi$ (rad)")
ax.set_ylabel(r"$\psi$ (rad)")
#ax.set_xlim(-np.pi/2,np.pi/2)
#ax.set_ylim(-np.pi/2,np.pi/2)
'''

# psi vs psi'

ax=fig.add_subplot(111)
ax.plot(datos[:,2], datos[:,3], linewidth=0.5, color="green")
ax.set_xlabel(r"$\psi$ (rad)")
ax.set_ylabel(r"$\dot{\psi}$ (rad)")


#psi vs phi'
"""
ax=fig.add_subplot(111)
ax.plot(datos[:,2], datos[:,1], linewidth=0.5, color="red")
ax.set_xlabel(r"$\psi$ (rad)")
ax.set_ylabel(r"$\dot{\psi}$ (rad)")
"""

fig.savefig(file_out)