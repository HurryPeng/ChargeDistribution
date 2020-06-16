import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

pointsFile = open("pointsPy.txt", mode = 'r')
points = eval(pointsFile.read())

x, y, z = [], [], []
for point in points:
    x.append(point[0])
    y.append(point[1])
    z.append(point[2])

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(x, y, z, s = 0.1)

plt.show()
