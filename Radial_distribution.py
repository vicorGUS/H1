import numpy as np
import matplotlib.pyplot as plt

pos = np.genfromtxt('pos_700.csv', delimiter=',')
pos_super = np.zeros([256*27, 3])

a = 4.13

"""for i in range(27):
    pos_super[256 * i:256 * (i + 1)][0] = pos[:][0] + a * np.floor(i/3)
    pos_super[256 * i:256 * (i + 1)][1] = pos[:][1] + a * np.floor((i+1)/3)
    pos_super[256 * i:256 * (i + 1)][2] = pos[:][2] + a * np.floor((i+2)/3)"""

dr = 0.2

#super_pos = np.stack([pos, np.array([pos[:, 0]+a, pos[:, 1], pos[:, 2]])])

pos100 = np.transpose(np.array([pos[:, 0]+a, pos[:, 1], pos[:, 2]]))
pos010 = np.transpose(np.array([pos[:, 0], pos[:, 1]+a, pos[:, 2]]))
pos001 = np.transpose(np.array([pos[:, 0], pos[:, 1], pos[:, 2]+a]))

pos200 = np.transpose(np.array([pos[:, 0]-a, pos[:, 1], pos[:, 2]]))
pos020 = np.transpose(np.array([pos[:, 0], pos[:, 1]-a, pos[:, 2]]))
pos002 = np.transpose(np.array([pos[:, 0], pos[:, 1], pos[:, 2]-a]))

pos110 = np.transpose(np.array([pos[:, 0]+a, pos[:, 1]+a, pos[:, 2]]))
pos011 = np.transpose(np.array([pos[:, 0], pos[:, 1]+a, pos[:, 2]+a]))
pos101 = np.transpose(np.array([pos[:, 0]+a, pos[:, 1], pos[:, 2]+a]))

pos220 = np.transpose(np.array([pos[:, 0]-a, pos[:, 1]-a, pos[:, 2]]))
pos022 = np.transpose(np.array([pos[:, 0], pos[:, 1]-a, pos[:, 2]-a]))
pos202 = np.transpose(np.array([pos[:, 0]-a, pos[:, 1], pos[:, 2]-a]))

pos102 = np.transpose(np.array([pos[:, 0]+a, pos[:, 1], pos[:, 2]-a]))
pos021 = np.transpose(np.array([pos[:, 0], pos[:, 1]-a, pos[:, 2]+a]))
pos210 = np.transpose(np.array([pos[:, 0]-a, pos[:, 1]+a, pos[:, 2]]))

pos201 = np.transpose(np.array([pos[:, 0]-a, pos[:, 1], pos[:, 2]+a]))
pos012 = np.transpose(np.array([pos[:, 0], pos[:, 1]+a, pos[:, 2]-a]))
pos120 = np.transpose(np.array([pos[:, 0]+a, pos[:, 1]-a, pos[:, 2]]))

pos112 = np.transpose(np.array([pos[:, 0]+a, pos[:, 1]+a, pos[:, 2]-a]))
pos121 = np.transpose(np.array([pos[:, 0]+a, pos[:, 1]-a, pos[:, 2]+a]))
pos211 = np.transpose(np.array([pos[:, 0]-a, pos[:, 1]+a, pos[:, 2]+a]))

pos221 = np.transpose(np.array([pos[:, 0]-a, pos[:, 1]-a, pos[:, 2]+a]))
pos212 = np.transpose(np.array([pos[:, 0]-a, pos[:, 1]+a, pos[:, 2]-a]))
pos122 = np.transpose(np.array([pos[:, 0]+a, pos[:, 1]-a, pos[:, 2]-a]))

pos111 = np.transpose(np.array([pos[:, 0]+a, pos[:, 1]+a, pos[:, 2]+a]))
pos222 = np.transpose(np.array([pos[:, 0]-a, pos[:, 1]-a, pos[:, 2]-a]))

super_pos = np.row_stack([pos, pos100, pos010, pos010, pos200, pos020, pos002, pos110, pos011, pos101, pos220, pos022,
                             pos202, pos102, pos021, pos210, pos201, pos012, pos120, pos112, pos121, pos211, pos221, pos212,
                             pos122, pos111, pos222])
print(super_pos.shape)

"""fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')

ax.scatter(super_pos[:,0], super_pos[:,1], super_pos[:,2])
plt.show()"""

deltas = super_pos[:, np.newaxis] - super_pos

kmax = 100
N = np.zeros([256, kmax])
N_ideal = np.zeros(kmax)

L = a * 4.0
V = L ** 3

"""for i in range(len(deltas)):
    for j in range(len(deltas)):
        for d in range(3):
            if deltas[i, j, d] > L/2:
                deltas[i, j, d] -= L
            if deltas[i, j, d] < -L/2:
                deltas[i, j, d] += L"""

dist = np.linalg.norm(deltas, axis=-1)

for i in range(256):
    for j in range(6912):
        for k in range(kmax):
            N_ideal[k] = 255 * 4 * np.pi / (3 * V) * (3 * k ** 2 - 3 * k + 1) * dr ** 3
            if (k - 1 / 2) * dr <= abs(dist[i, j]) < k * dr:
                N[i][k] += 1

N_avg = np.zeros(kmax)

for i in range(256):
    for k in range(kmax):
        N_avg[k] += 1/256 * N[i][k]

g = N_avg / N_ideal
k_line = np.linspace(0, 10, kmax)

plt.plot(k_line, N_avg / N_ideal)

plt.show()
"""
plt.plot(k_line, N_ideal)
plt.show()"""


