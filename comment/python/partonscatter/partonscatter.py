#!/usr/bin/env python
# _*_ coding: utf-8 _*_

"""document of the module"""


import os
import sys
import subprocess
import matplotlib.pyplot as plt
import numpy as np


fileobj = open("/storage/third/dataForAmpt/final/phase/cluster4/ampt/testOut/partonxiaohai.dat")
alllines = fileobj.readlines()

colorList = ["b", "g", "r", "c", "m", "y", "k"]
markerList = [",", "o", "v", "^", "<", ">", "*", "+"]
markerBig = [20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 60, 75]
countColor = 0
fig = plt.figure(figsize=(10, 10), dpi=100)
xlist = []
ylist = []
zlist = []

for line in alllines:
    if len(line.strip("\n").split()) == 1:
        xlist = []
        ylist = []
        zlist = []
        countColor += 1
        continue
    xlist.append(float(line.strip("\n").split()[0]))
    ylist.append(float(line.strip("\n").split()[1]))
    zlist.append(float(line.strip("\n").split()[2]))
    colorstyle = colorList[countColor%7]
    markerStyle = markerList[countColor%8]
    bigStyle = markerBig[countColor%12]
    plt.scatter(xlist, ylist, s=bigStyle, c=colorstyle, marker=markerStyle)

plt.xlim(-6, 6)
plt.ylim(-6, 6)
plt.savefig("parton.png")
plt.show()
