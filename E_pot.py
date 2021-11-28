import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('E_pot.csv', delimiter=',', skip_header=1)

plt.figure(figsize=(16, 10))

plt.scatter(array[:, 0], array[:, 1], s=300)
plt.plot(array[:, 0], array[:, 1], '--')

f = np.poly1d(np.polyfit(array[:, 0], array[:, 1], deg=2))
fline = np.linspace(63.5, 68.0, 100)
a0 = np.argmin(f(fline))

print("Teoretical lattice parameter:", fline[a0]**(1/3))

plt.grid()
plt.xlabel(r'Volume (Å${^3}$)')
plt.ylabel('Energy (eV/unit cell)')

plt.savefig('E_pot.jpg')