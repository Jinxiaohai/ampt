#!/usr/bin/env python
# _*_ coding: utf-8 _*_

"""document of the module"""


import os
import sys
import subprocess
import matplotlib.pyplot as plt
import numpy as np


fileObj = open("../ana/zpc.dat")
d, u, s, dbar, ubar, sbar = 0, 0, 0, 0, 0, 0
event = 0
while True:
    headLine = fileObj.readline()
    if not headLine:
        break
    event += 1
    multi = int(headLine.strip().split()[2])
    for eachTrack in range(multi):
        id = int(fileObj.readline().split()[0])
        if id == 1:
            d += 1
        if id == 2:
            u += 1
        if id == 3:
            s += 1
        if id == -1:
            dbar += 1
        if id == -2:
            ubar += 1
        if id == -3:
            sbar += 1

print "d is ", d, event, float(d)/event
print "u is ", u, event, float(u)/event
print "s is ", s, event, float(s)/event
print "dbar is", dbar, event, float(dbar)/event
print "ubar is", ubar, event, float(ubar)/event
print "sbar is", sbar, event, float(sbar)/event
