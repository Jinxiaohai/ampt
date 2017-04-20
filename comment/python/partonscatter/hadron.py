#!/usr/bin/env python
# _*_ coding: utf-8 _*_

"""document of the module"""


import os
import sys
import subprocess
import matplotlib.pyplot as plt
import numpy as np


fileobj = open("../../testOut/hadroncoordinatexiaohai.dat")
alllines = fileobj.readlines()

xlist = []
ylist = []
zlist = []

for line in alllines:
    xlist.append(float(line.strip("\n").split()[0]))
    ylist.append(float(line.strip("\n").split()[1]))
    zlist.append(float(line.strip("\n").split()[2]))
    
fig = plt.figure(figsize=(10, 10), dpi=100)
plt.scatter(xlist, ylist, s=40, c='g', marker='o')
plt.xlim(-10, 10)
plt.ylim(-10, 10)
plt.savefig("hadron.png")
plt.show()
