from matplotlib import pyplot as plt
import numpy as np

file_in = "lyapunov4.dat"
file_out = "lyapunov4"

datos=np.loadtxt(file_in)

fig=plt.figure()


ax=fig.add_subplot(111)
ax.errorbar(np.log10(datos[:,0]),datos[:,1], yerr=datos[:,2], color="green", ecolor="red", fmt="o-", capsize=2, markersize=1.5, linewidth=0.5, elinewidth=0.3, label=r'$\lambda_\phi$')
ax.errorbar(np.log10(datos[:,0]),datos[:,3], yerr=datos[:,4], color="blue", ecolor="red", fmt="o-", capsize=2, markersize=1.5, linewidth=0.5, elinewidth=0.3, label=r'$\lambda_\psi$')
ax.set_xlabel(r'$log_{10}(d_0)$')
ax.set_ylabel("Coeficiente de Lyapunov")
ax.legend()

fig.savefig(file_out)
