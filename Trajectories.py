import numpy as np
import matplotlib.pyplot as plt

trajectories = np.genfromtxt('q.csv', delimiter=',')

plt.figure(figsize=(16, 10))

plt.plot(trajectories[:, 0], trajectories[:, 1], label='Particle 1')
plt.plot(trajectories[:, 2], trajectories[:, 3], label='Particle 2')
plt.plot(trajectories[:, 4], trajectories[:, 5], label='Particle 3')

plt.legend(loc='best')
plt.xlabel("Position (Å)")
plt.ylabel("Position (Å)")

plt.savefig('Trajectories.jpg')