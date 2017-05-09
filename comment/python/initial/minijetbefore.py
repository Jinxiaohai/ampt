#!/usr/bin/env python
# _*_ coding: utf-8 _*_

"""document of the module"""


import os
import sys
import subprocess
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, Circle


fileobj = open("./minijet-initial-beforePropagation.dat")
xlist = []
ylist = []
jet = fileobj.readlines()


projectileX = []
projectileY = []
targetX = []
targetY = []
njsgsX = []
njsgsY = []

for inum in jet:
    flag = int(inum.split()[9])
    x = float(inum.split()[5])
    y = float(inum.split()[6])
    if flag == 1:
        projectileX.append(x)
        projectileY.append(y)
    if flag == 2:
        targetX.append(x)
        targetY.append(y)
    if flag == 3:
        njsgsX.append(x)
        njsgsY.append(y)


fig = plt.figure(figsize=(10,10), dpi=150)
fig.set_size_inches(6, 6)

plt.scatter(njsgsX, njsgsY, s=100, c='w', marker='o')
ax1 = plt.scatter(projectileX, projectileY, s=100, c='r', marker='o')
ax2 = plt.scatter(targetX, targetY, s=100, c='g', marker='o')

plt.xlim(-10, 12)
plt.ylim(-10, 10)
plt.xlabel("x [fm]")
plt.ylabel("y [fm]")
plt.legend((ax1, ax2), ("Au", "Cu"))
plt.show()
