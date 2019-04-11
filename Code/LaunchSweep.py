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
    sys.argv[1:]=["3",]   # Default for NprocessMax = 1
    
NprocessMax=int(sys.argv[1])

WaitTime = 10


Data=load("SweepLists/SweepList.npz")
SeriesList=Data["SeriesList"] 



L=int(SeriesList[0][0][0])
Nperiods = int(SeriesList[0][0][1])
PBC=int(SeriesList[0][0][2])

NSeries = len(SeriesList)

Nprocess = 0
plist=[]
termlist=[]
def PrintStatus():
                                

    print(f"{len(plist)-len(termlist)} processes running, {len(termlist)} terminated")
    print("")
    print(f"    Running series     : {list(where([q.is_alive() for q in plist])[0])}")
    print(f"    Terminated series  : {list(where([not q.is_alive() for q in plist])[0])}")


def CheckActiveProcesses(plist):
    Nalive = 0
    
    for p in plist:
        if p.is_alive():
            Nalive+=1
        else:
            # Join dead processes to avoid zombies on server. 
            if not p in termlist:
                p.join()
                termlist.append(p)                
                print("")
                print(f"========================= Series {plist.index(p)} terminated ========================")                 
                print("")
                PrintStatus()
                print("")
                print("----------------------------------------------------------------------")
                print("")      
                
    sys.stdout.flush()
    sys.stderr.flush()

    return Nalive 

for ns in range(0,NSeries):
    

    Series=SeriesList[ns]
    
    # Sweep file saves parameters as strings, to avoid rounding errors. 
    Series=[[int(x[0]),int(x[1]),int(x[2]),float(x[3]),float(x[4])] for x in Series]
    
    time.sleep(WaitTime/100)
    
    FileString=f"{int(L+0.1)}_{int(Nperiods+0.1)}_{int(PBC+0.1)}"
    
    
    while True:
        
        ActiveProcesses=CheckActiveProcesses(plist)
#        print("")
#        print("----------------------------------------------------------------------")
#        print(f"{ActiveProcesses} processes running, {len(plist)-ActiveProcesses} terminated")
#        print("----------------------------------------------------------------------")
#        print("")
        if ActiveProcesses<NprocessMax:
            break
        else:
            sys.stdout.flush()
            sys.stderr.flush()

            time.sleep(WaitTime)
        

    print("")
    print(f"========================= Launching series {ns} =========================")
    print("")


    p = Process(target = PF.EnergyAbsorptionSweep,args=(Series,FileString),kwargs={"OP":True,"PreString":f"     Series {ns}: "})
   
    p.start()

    plist.append(p)
    


while True:
    
    ActiveProcesses=CheckActiveProcesses(plist)
    

    if ActiveProcesses==0:
        break
    else:

        time.sleep(WaitTime)

sys.stdout.flush()
sys.stderr.flush()
print("")
print("Done.")








