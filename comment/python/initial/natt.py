#!/usr/bin/env python
# _*_ coding: utf-8 _*_

"""document of the module"""


import os
import sys
import subprocess
import matplotlib.pyplot as plt
import numpy as np


file = open("../../testOut/natt2xiaohai.dat")

nsgx = []
nsgy = []
npjx = []
npjy = []
ntjx = []
ntjy = []

for eachline in file:
    istrg = int(eachline.split()[9])
    x = float(eachline.split()[5])
    y = float(eachline.split()[6])
    if istrg > 0 and istrg < 10000:
        nsgx.append(x)
        nsgy.append(y)
    if istrg > 10000 and istrg < 20000:
        npjx.append(x)
        npjy.append(y)
    if istrg > 20000:
        ntjx.append(x)
        ntjy.append(y)

fig = plt.figure(figsize=(10,10), dpi=100)
plt.scatter(nsgx, nsgy, s=40, marker='o', c='w')
plt.scatter(npjx, npjy, s=40, marker='o', c='r')
plt.scatter(ntjx, ntjy, s=40, marker='o', c='g')
plt.xlim(-10, 10)
plt.ylim(-10, 10)
plt.xlabel("x [fm]")
plt.ylabel("y [fm]")
plt.show()
