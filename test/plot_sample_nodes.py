import numpy as np
import matplotlib.pyplot as plt
from itertools import product, combinations

from os import path

CSD = path.dirname(path.abspath(__file__))

sample_nodes = np.load(path.join(CSD, 'sample_nodes.npy'))

# Plot a cube with unit length
# draw a cube with unit length (1,1,1)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
r = [0, 1]
for s, e in combinations(np.array(list(product(r, r, r))), 2):
   if np.sum(np.abs(s-e)) == r[1]-r[0]:
      ax.plot3D(*zip(s, e), color="green")


ax.legend()
sf = 20.0
c = ['b'] * 8 + ['g'] * 24 + ['k'] * 18
s = [2.0 * sf] * 8 + [1.0 * sf] * 24 + [0.5 * sf] * 18
ax.scatter(sample_nodes[:, 0], sample_nodes[:, 1], sample_nodes[:, 2], c=c, s=s)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()
