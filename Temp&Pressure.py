import numpy as np
import matplotlib.pyplot as plt

array_T = np.genfromtxt('T_700.csv', delimiter=',', skip_header=1)
array_P = np.genfromtxt('P_700.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots(1, 2)

ax[0].plot(array_T[20:, 0], array_T[20:, 1], color='r')
ax[0].set_ylabel('Temperature [K]')

ax[1].plot(array_P[20:, 0], array_P[20:, 1], color='b')
ax[1].set_ylabel('Pressure [MPa]')

plt.tight_layout()

fig.savefig('Temp&Pressure_700.jpg')

plt.figure(figsize=(16,10))

plt.plot(array_T[20:, 0], array_T[20:, 1], color='r')
plt.ylabel('Temperature [K]')
plt.xlabel('Time (ps)')

plt.savefig('Temperature_700.jpg')

plt.figure(figsize=(16,10))

plt.plot(array_P[20:, 0], array_P[20:, 1], color='b')
plt.ylabel('Pressure [MPa]')
plt.xlabel('Time (ps)')

plt.savefig('Pressure_700.jpg')