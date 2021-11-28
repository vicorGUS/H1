import numpy as np
import matplotlib.pyplot as plt

array_p = np.genfromtxt('Ep_700.csv', delimiter=',', skip_header=1)
array_k = np.genfromtxt('Ek_700.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots()

ax.plot(array_p[:, 0], array_p[:, 1] / 64, label='Potential Energy')
ax.plot(array_k[:, 0], array_k[:, 1] / 64, label='Kinetic Energy')
ax.plot(array_k[:, 0], (array_k[:, 1] + array_p[:, 1]) / 64, label='Total Energy')

ax.set_xlabel('Time (ps)')
ax.set_ylabel('Energy (eV/unit cell)')
ax.grid()

plt.legend(loc='best')

fig.savefig('E_700.jpg')