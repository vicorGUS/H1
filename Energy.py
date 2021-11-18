import numpy as np
import matplotlib.pyplot as plt

#array = np.genfromtxt('E_pot.csv', delimiter=',', skip_header=1)
array_p = np.genfromtxt('Ep.csv', delimiter=',', skip_header=1)
array_k = np.genfromtxt('Ek.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots()

#ax.plot(array[:, 0], array[:, 1], '--')
ax.plot(array_p[1:, 0], array_p[1:, 1], label='Potential Energy')
ax.plot(array_k[1:, 0], array_k[1:, 1], label='Kinetic Energy')
ax.plot(array_k[1:, 0], array_k[1:, 1] + array_p[1:, 1], label='Total Energy')
#ax.scatter(array[:, 0], array[:, 1])

ax.set_xlabel('Time (ps)')
ax.set_ylabel('Energy (eV/unit cell)')
ax.grid()
plt.legend(loc='best')

fig.savefig('E.jpg')