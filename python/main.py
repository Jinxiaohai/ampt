#!/usr/bin/env python
# _*_ coding: utf-8 _*_

"""document of the module"""


import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import re


def CompareStack(fileline):
    """
    return the list about the subroutine
    """
    countModuleNum = 0
    stack = []
    content = []
    allList = []
    for eachline in fileline:
        sourceFileLine = eachline             # delete the '\n'
        eachline = eachline.strip()
        # print eachline
        if ((eachline.split()[0].lower() == "subroutine")
            or (eachline.split()[0].lower() == "program")
            or (eachline.split()[0].lower() == "function")
            or (eachline.split()[0].lower() == "real"
                and eachline.split()[1].lower() == "function")
            or (eachline.split()[0].lower() == "integer"
                and eachline.split()[1].lower() == "function")
            or (eachline.split()[0].lower() == "block"
                and eachline.split()[1].lower() == "data")
            or (eachline.split()[0].lower() == "double"
                and eachline.split()[1].lower() == "precision"
                and eachline.split()[2].lower() == "function")):
            stack.append("start")
            countModuleNum += 1
            content = []
        content.append(sourceFileLine)
        if eachline.lower() == "end":
            stack.pop()
            if len(stack) % 2 == 0:
                # print 'OK OK OK !!!'
                # print
                allList.append(content)
    return allList


def GetName(eachSubroutine):
    if eachSubroutine[0].split()[0].lower() == "subroutine":
        fileName = eachSubroutine[0].split()[1].lower().split('(')[0]
    elif eachSubroutine[0].split()[0].lower() == "program":
        fileName = eachSubroutine[0].split()[1].lower().split('(')[0]
    elif (eachSubroutine[0].split()[0].lower() == "function"):
        fileName = eachSubroutine[0].split()[1].lower().split('(')[0]
    elif (eachSubroutine[0].split()[0].lower() == "real"
          and eachSubroutine[0].split()[1].lower() == "function"):
        fileName = eachSubroutine[0].split()[2].lower().split('(')[0]
    elif (eachSubroutine[0].split()[0].lower() == "integer"
          and eachSubroutine[0].split()[1].lower() == "function"):
        fileName = eachSubroutine[0].split()[2].lower().split('(')[0]
    elif (eachSubroutine[0].split()[0].lower() == "block"
          and eachSubroutine[0].split()[1].lower() == "data"):
        fileName = eachSubroutine[0].split()[2].lower().split('(')[0]
    elif (eachSubroutine[0].split()[0].lower() == "double"
          and eachSubroutine[0].split()[1].lower() == "precision"
          and eachSubroutine[0].split()[2].lower() == "function"):
        fileName = eachSubroutine[0].split()[3].lower().split('(')[0]
    return fileName


def dealWithList(alllist):
    for eachSubroutine in alllist:
        fileName = ""
        if eachSubroutine[0].split()[0].lower() == "subroutine":
            fileName = eachSubroutine[0].split()[1].lower().split('(')[0]
        elif eachSubroutine[0].split()[0].lower() == "program":
            fileName = eachSubroutine[0].split()[1].lower().split('(')[0]
        elif (eachSubroutine[0].split()[0].lower() == "function"):
            fileName = eachSubroutine[0].split()[1].lower().split('(')[0]
        elif (eachSubroutine[0].split()[0].lower() == "real"
              and eachSubroutine[0].split()[1].lower() == "function"):
            fileName = eachSubroutine[0].split()[2].lower().split('(')[0]
        elif (eachSubroutine[0].split()[0].lower() == "integer"
              and eachSubroutine[0].split()[1].lower() == "function"):
            fileName = eachSubroutine[0].split()[2].lower().split('(')[0]
        elif (eachSubroutine[0].split()[0].lower() == "block"
              and eachSubroutine[0].split()[1].lower() == "data"):
            fileName = eachSubroutine[0].split()[2].lower().split('(')[0]
        elif (eachSubroutine[0].split()[0].lower() == "double"
              and eachSubroutine[0].split()[1].lower() == "precision"
              and eachSubroutine[0].split()[2].lower() == "function"):
            fileName = eachSubroutine[0].split()[3].lower().split('(')[0]
        # try:
        #     writeName = open("./amptSource/" + fileName + ".f", 'w')
        # except IOError:
        #     pass
        # for eachline in eachSubroutine:
        #     writeName.writelines(eachline + '\n')
        # writeName.close()


def dealWith(allList):
    ROOT = []
    print GetName(allList[0])
    ROOT.append(signleSubroutine(allList[0], allList))
    return ROOT


def signleSubroutine(subroutine, allList):
    flag = False
    theName = GetName(subroutine)
    saveName = []
    saveName.append(theName)
    for eachline in subroutine:
        eachlineList = eachline.split()
        for i in range(len(eachlineList)):
            if eachlineList[i].lower() == "call":
                flag = True
                continue
            if flag is True:
                subroutineName = eachlineList[i].lower().split('(')[0]
                for j in allList:
                    if GetName(j) == subroutineName:
                        signleSubroutine(j, allList)
                flag = False
    return saveName


def GetEachline(fileobj):
    """
    fileobj(had deleted the comment line)!!!
    """
    fileline = []
    for eachline in fileobj:
        kongline = eachline.strip('\n').strip()
        if len(kongline) != 0:
            fileline.append(eachline.strip('\n'))
    return fileline


def GetSubroutine(fileline):
    subroutinelist = []
    for eachline in fileline:
        eachLineList = eachline.split()
        if eachLineList[0].lower() == 'subroutine':
            subroutinelist.append(eachLineList[1].split('(')[0])
    return subroutinelist


def GetCallSubroutine(fileline):
    tag = False
    for eachline in fileline:
        eachlineList = eachline.split()
        for i in range(len(eachlineList)):
            if eachlineList[i].lower() == 'call':
                tag = True
                continue
            if tag is True:
                print eachlineList[i].split("(")[0],
                tag = False


def main(file):
    fileobj = open(file)
    fileline = GetEachline(fileobj)
    # allList = CompareStack(fileline)
    # dealWithList(allList)
    # dealWith(allList)
    subroutine = GetSubroutine(fileline)
    # GetCallSubroutine(fileline)
    print
    for i in subroutine:
        print i
    fileobj.close()


if __name__ == '__main__':
    # filename = raw_input("input file name: ")
    main("./log.f")
