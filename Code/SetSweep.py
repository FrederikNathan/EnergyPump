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
import DataAnalysisFunctions as D
from multiprocessing import Process as Process 

from subprocess import call



    
    
SweepListName="SweepList.npz"

# =============================================================================
# 1: Set parameters
# =============================================================================

#L=400
#Nperiods =10000
#PBC=0
#
#Nseries=500
#
#OmegaMin=1*2*pi
#OmegaMax=3*2*pi
#
#Wmin=0*2*pi
#Wmax=3*2*pi
#
#dOmega = 0.1*1.6183
#dW=0.1
#
#Wrange = arange(Wmin+dW,Wmax,dW)
#OmegaRange = arange(OmegaMin+dOmega,OmegaMax,dOmega)
#
#
#NSamplesMin=2  # Minimal number of samples

# =============================================================================
# 2: Sweep list commands
# =============================================================================
""" Here goes the command generating the sweeplist

Examples:
    
For sequence of omega-sweeps, use

SeriesList=[]
for W in Wrange:
    
     Runlist=D.OmegaSweep(L,Nperiods,PBC,0,OmegaRange)
     SeriesList.append(Series)
    

For filling sweep, use

RunList=D.FillingSweep(L,Nperiods,PBC,Wrange,OmegaRange,NsamplesMin)
SeriesList = D.SeriesDivider(Runlist,Nseries)
"""
Wrange = array([0.5*2*pi])

dOmega=0.1*1.61823
OmegaRange=2*pi*arange(0.8,4,dOmega)

Lrange=[7,8,9]#[25,50,100,200,400,800]
NperiodsRange =[10,20,30]# [10,50,100,500,1000,5000,10000,50000]

PBC=0
Ncopies = 1
NprocessMax=10
# For sequence of omega-sweeps
RunList=[]
for Nperiods in NperiodsRange:
    
    for L in Lrange:
        
        for W in Wrange:
            
            for nc in range(0,Ncopies):
                RunList=RunList+D.OmegaSweep(L,Nperiods,PBC,W,OmegaRange)

RunList=[RunList[n] for n in permutation(arange(0,len(RunList)))]
SeriesList=D.SeriesDivider(RunList,NprocessMax)
Nruns = sum([len(x) for x in SeriesList])
##
##
#omega1=2*pi
#
#Wrange=arange(1.5*2*pi,2.2*2*pi,0.1)
#
#dOmega=0.1*1.61823
#OmegaRange=arange(0.0*omega1+dOmega,omega1*2.5,dOmega)
#
###LNPlist[10000,20000,40000,80000]
##Lrange=[400]
##NperiodsRange = [10000]
##
##PBC=0
##Ncopies = 2
##
### For sequence of omega-sweeps
##SeriesList=[]
##for Nperiods in NperiodsRange:
##    
##    for L in Lrange:
##        
##        for W in Wrange:
##            for nc in range(0,Ncopies):
##                Series=D.OmegaSweep(L,Nperiods,PBC,W,OmegaRange)
##                SeriesList.append(Series)
##                
##Nruns = sum([len(x) for x in SeriesList])
##
#NSamplesMin=1
#Nseries=50
#Nperiods = 10000
#L=400
#PBC=0
#
#
#Wrange=arange(1.5*2*pi,2.2*2*pi,0.1)
#
#dOmega=0.1*1.61823
#OmegaRange=arange(0.0*omega1+dOmega,omega1*2.5,dOmega)
#
## For filling sweep
#RunList=D.FillingSweep(L,Nperiods,PBC,Wrange,OmegaRange,NSamplesMin)
#SeriesList = D.SeriesDivider(RunList,Nseries)
#Nruns = len(RunList)
print("")
print(f"Number of runs   : {Nruns}   ({len(SeriesList)} series)")
print("")

# =============================================================================
# Divide into series and save
# =============================================================================

        
savez(f"SweepLists/{SweepListName}",SeriesList=SeriesList)
    
    

    
    
    