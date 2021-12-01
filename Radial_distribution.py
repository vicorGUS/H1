import numpy as np
import matplotlib.pyplot as plt

pos = np.genfromtxt('pos_700.csv', delimiter=',')
pos_super = np.zeros([256*27, 3])

for i in range(27):
    pos_super[256 * i:256 * (i + 1)][0] = pos[:][0] + 4.13 * np.floor(i/3)
    pos_super[256 * i:256 * (i + 1)][1] = pos[:][1] + 4.13 * np.floor((i+1)/3)
    pos_super[256 * i:256 * (i + 1)][2] = pos[:][2] + 4.13 * np.floor((i+2)/3)

print(pos_super.shape)

dr = 0.1

deltas = pos_super[:, np.newaxis] - pos_super

kmax = 300
N = np.zeros([256, kmax])
N_ideal = np.zeros(kmax)

L = 4.13 * 4.0
V = L ** 3

"""for i in range(len(deltas)):
    for j in range(len(deltas)):
        for d in range(3):
            if deltas[i, j, d] > L/2:
                deltas[i, j, d] -= L
            if deltas[i, j, d] < -L/2:
                deltas[i, j, d] += L"""

dist = np.linalg.norm(deltas, axis=-1)

avg = np.zeros(kmax)

for i in range(256):
    for j in range(6912):
        for k in range(kmax):
            N_ideal[k] = 255 * 4 * np.pi / (3 * V) * (3 * k ** 2 - 3 * k + 1) * dr ** 3
            if (k - 1 / 2) * dr <= abs(dist[i+3327, j]) < k * dr:
                N[i][k] += 1

N_avg = np.zeros(kmax)

for i in range(256):
    for k in range(kmax):
        N_avg[k] += 1/256 * N[i][k]

g = N_avg / N_ideal
k_line = np.linspace(0, 10, kmax)

plt.plot(k_line, N_avg)

plt.show()

"""plt.plot(k_line, N_ideal)
plt.show()"""


