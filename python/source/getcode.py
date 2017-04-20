#!/usr/bin/env python
# _*_ coding: utf-8 _*_

"""document of the module"""


import re


fileobj = open("./test.f")

fileline = []

for eachline in fileobj:
    fileline.append(eachline.strip('\n'))


for eachline in fileline:
    p = re.compile(r'[^c,^*,^C]')
    result = p.match(eachline)
    if result:
        print eachline
    
