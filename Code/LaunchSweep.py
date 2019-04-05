#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 14:39:24 2019

@author: frederik
"""

# =============================================================================
# 1: Parameters
# =============================================================================

import sys
import os 
import os.path as op

from scipy import *
from scipy.linalg import *
from numpy import ufunc as ufunc
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numpy.random import *
from scipy.misc import factorial as factorial

import datetime
import logging as l
import time

import TightBindingModel as TBM 

import LatticeModel as LM

import BasicFunctions as BF

import PumpFunctions as PF

from multiprocessing import Process as Process 

if len(sys.argv)==1:
    sys.argv[1:]=["1",]   # Default for NprocessMax = 1
    print("Hej")
    
NprocessMax=int(sys.argv[1])

WaitTime = 1


Data=load("SweepLists/SweepList.npz")
SeriesList=Data["SeriesList"] 



L=SeriesList[0][0][0]
Nperiods = SeriesList[0][0][1]
PBC=SeriesList[0][0][2]

NSeries = len(SeriesList)

Nprocess = 0
plist=[]



def CheckActiveProcesses(plist):
    Nalive = 0
    
    for p in plist:
        if p.is_alive():
            Nalive+=1
    

    return Nalive 

for ns in range(0,NSeries):
    

    Series=SeriesList[ns]
    time.sleep(0.2)
    
    FileString=f"{int(L+0.1)}_{int(Nperiods+0.1)}_{int(PBC+0.1)}"
    
    
    while True:
        
        ActiveProcesses=CheckActiveProcesses(plist)
        print(f"    {ActiveProcesses} processes running, {len(plist)-ActiveProcesses} terminated")
        print("")
        if ActiveProcesses<=NprocessMax:
            break
        else:

            time.sleep(WaitTime)
        


        
    p = Process(target = PF.EnergyAbsorptionSweep,args=(Series,FileString))
   
    p.start()

    print("Launching 1 process")
    print("")
    plist.append(p)














