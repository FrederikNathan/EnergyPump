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
Wrange = 2*pi*array([0.3,0.8,1.2,2.3])

dOmega=0.1*1.61823
OmegaRange=2*pi*arange(0,7,dOmega)

Lrange=[20,50,100,200,400,800]
NperiodsRange = [10,100,1000,10000,50000]

PBC=0
Ncopies = 5

# For sequence of omega-sweeps
SeriesList=[]
for Nperiods in NperiodsRange:
    
    for L in Lrange:
        
        for W in Wrange:
            
            for nc in range(0,Ncopies):
                Series=D.OmegaSweep(L,Nperiods,PBC,W,OmegaRange)
                SeriesList.append(Series)
                
Nruns = sum([len(x) for x in SeriesList])

## For filling sweep
#RunList=D.FillingSweep(L,Nperiods,PBC,Wrange,OmegaRange,NSamplesMin)
#SeriesList = D.SeriesDivider(RunList,Nseries)
#Nruns = len(RunList)

print(f"Number of runs : {Nruns}")
print("")

# =============================================================================
# Divide into series and save
# =============================================================================

        
savez(f"Sweeplists/{SweepListName}",SeriesList=SeriesList)
    
    

    
    
    