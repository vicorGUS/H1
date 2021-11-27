import numpy as np
import matplotlib.pyplot as plt

array_T = np.genfromtxt('T.csv', delimiter=',', skip_header=1)
array_P = np.genfromtxt('P.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots(1, 2)

ax[0].plot(array_T[:, 0], array_T[:, 1], color='r')
#ax[0].set_xticks((0.0, 0.001, 0.002, 0.003, 0.004, 0.005))
#ax[0].set_xticklabels(('0.0', r'$1\tau_T$', r'$2\tau_T$', r'$3\tau_T$', r'$4\tau_T$', r'$5\tau_T$'))
ax[0].set_ylabel('Temperature [K]')

ax[1].plot(array_P[:, 0], array_P[:, 1], color='b')
#ax[1].set_xticks((0.0250, 0.0275, 0.0300, 0.0325, 0.0350))
#ax[1].set_xticklabels(('0.0', r'$1\tau_P$', r'$2\tau_P$', r'$3\tau_P$', r'$4\tau_P$'))
ax[1].set_ylabel('Pressure [MPa]')

plt.tight_layout()

fig.savefig('Temp&Pressure.jpg')

plt.figure(figsize=(16,10))

plt.plot(array_T[:, 0], array_T[:, 1], color='r')
plt.ylabel('Temperature [K]')
plt.xlabel('Time (ps)')

plt.savefig('Temperature.jpg')

plt.figure(figsize=(16,10))

plt.plot(array_P[:, 0], array_P[:, 1], color='b')
plt.ylabel('Pressure [MPa]')
plt.xlabel('Time (ps)')

plt.savefig('Pressure.jpg')