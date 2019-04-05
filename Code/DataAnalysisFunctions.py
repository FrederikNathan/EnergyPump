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

DataDir = "../Data/"
CleanDataDir="../DataClean/"

FileList = os.listdir(DataDir)

def IsNPZ(FileName):
    if len(FileName)>4:
        if FileName[-4:]==".npz":
            return True
        else:
            return False
    else:
        return False
    

def InputInterpreter(Word):
    try:
        Out = int(Word)
    except ValueError:
        Out = Word
            
    return Out
def FileNameParser(FileName):
    
    OutList=[]
    Word = ""
    for s in FileName:
        if not s == "_":
            Word+=s
            
        else:

            OutList.append(InputInterpreter(Word))
            Word=""
            
    OutList.append(InputInterpreter(Word))

    
    return OutList       
            
        
def FindIndex(Value,Range):
    dR=Range[1]-Range[0]
    y=argmin(abs(Range-Value))
    if amin(abs(Range-Value))>dR/2:
        return None
    else:
        return y 
    
    
# =============================================================================
# Sweep generators 
# =============================================================================
def FillingSweep(L,Nperiods,PBC,Wrange,OmegaRange,NSamplesMin):
    # =============================================================================
    # 2: Set and define variables 
    # =============================================================================
    
    
    NW=len(Wrange)
    NO=len(OmegaRange)
    
    Count = zeros((NW,NO),dtype=int)
    
    FileName = f"{L}_{Nperiods}_{PBC}.npz"
    Path = "../DataClean/"+FileName
    
    
    
    
    
    
    # =============================================================================
    # 3: Determine which parameters to run with
    # =============================================================================
    
    
    if os.path.isfile(Path):
        Data=load(Path)
        Wlist=Data["Wlist"]
        Om2list=Data["Om2list"]
        Ndata = len(Wlist)
        Elist=Data["Elist"]
        
        
        
        for n in range(0,Ndata):
            w = Wlist[n]
            o = Om2list[n]
            
            nOm=FindIndex(o,OmegaRange)
            nW=FindIndex(w,Wrange)
            
            if nOm!=None and nW!= None:
                
                Count[nW,nOm]+=1
        
    
    
        
    Runs = NSamplesMin-Count
    Runs=Runs*(Runs>0).astype(int)
    
    Runlist=[]
    for nw in range(0,NW):
        for no in range(0,NO):
            W=Wrange[nw]
            Om = OmegaRange[no]
            
            RunParameters =[L,Nperiods,PBC,Om,W]
            
            NRuns = Runs[nw,no]
            Runlist=Runlist+[RunParameters for x in range(0,NRuns)]
            
    
    

    
    return Runlist




def WSweep(L,Nperiods,PBC,Wrange,Omega2):
    Runlist=[]
    for W in Wrange:
        
        RunParameters =[L,Nperiods,PBC,Omega2,W]
        
        Runlist.append(RunParameters)
                
    
    return Runlist

def OmegaSweep(L,Nperiods,PBC,W,Omega2Range):
    Runlist=[]
    for Omega2 in Omega2Range:
        
        RunParameters =[L,Nperiods,PBC,Omega2,W]
        
        Runlist.append(RunParameters) # for x in range(0,NRuns)]
                
    
    return Runlist

def SeriesDivider(RunList,NprocessMax,MaxSeriesLength=10^4):
    """ Divide runlist into series for paralellization"""
    Nruns = len(RunList)
    SeriesSize = min(MaxSeriesLength,int(Nruns/NprocessMax/2)+1) # Distribute runs evenly    
    NSeries = int(Nruns /SeriesSize)+1
    
    OutList=[]
    

    for ns in range(0,NSeries):
        
        n0= min(ns*SeriesSize,Nruns)
        n1 = min((ns+1)*SeriesSize,Nruns)
    
        Series=[RunList[n] for n in range(n0,n1)]
        
    
        
        OutList.append(Series)
    
    return OutList
        