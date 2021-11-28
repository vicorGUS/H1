import numpy as np

k = 8.6173e-5
N = 256
T_500 = 500 + 275.15
T_700 = 700 + 275.15

Ek_500 = np.genfromtxt('Ek_500.csv', delimiter=',', skip_header=1)
Ek_700 = np.genfromtxt('Ek_700.csv', delimiter=',', skip_header=1)

dEk_500 = np.mean(Ek_500[15000:, 1]**2) - np.mean(Ek_500[15000:, 1])**2
dEk_700 = np.mean(Ek_700[15000:, 1]**2) - np.mean(Ek_700[15000:, 1])**2

Cv_500 = (3 * N * k / 2) * (1 - (2 / (3 * N * k**2 * T_500**2)) * dEk_500) ** -1
Cv_700 = (3 * N * k / 2) * (1 - (2 / (3 * N * k**2 * T_700**2)) * dEk_700) ** -1

print(Cv_500)
print(Cv_700)
