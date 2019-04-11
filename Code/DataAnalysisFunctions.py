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
    SeriesSize = min(MaxSeriesLength,int(Nruns/NprocessMax)+1) # Distribute runs evenly    
    NSeries = int(Nruns /SeriesSize)+1
    
    OutList=[]
    

    for ns in range(0,NSeries):
        
        n0= min(ns*SeriesSize,Nruns)
        n1 = min((ns+1)*SeriesSize,Nruns)
    
        Series=[RunList[n] for n in range(n0,n1)]
        
    
        
        OutList.append(Series)
    
    return OutList

# =============================================================================
# Data extractors
# =============================================================================
def FindData(Criterion):
    """ Find data that match the crieterion. 
    
    Criterion: callable, of the form f(Nperiods,L,pbc,Wlist,OmList). 
    Should return Indices (array of ints), of entries in Wlist,Omlist, where criterion is met. 
    Should return empty array (array([])) if no data points match the criterion. 
    
    returns
    
    [Earray,CmeanArray,OmArray,Warray,LegendList]
    
    Earray: ordered array with energy abpsorptions of the data that match the criterion, divided after the values of L,Nperiods and PBC
    Cmeanarray: mean correlation length
    OmArray: frequencies of mode 2
    Warray: disorder strengths
    SizeList: size of systems. Earray[n] is a list of data taken with [L,Nperiods,PBC]=SizeList[n]
    
    The data are sorted such that 
    Earray[n][z], CmeanArray[n][z] are obtained with parameters [L,Nperiods,PBC]=SizeList[n], W=Warray[n][z], Om2=Om2Array[n][z]
    """
    FileList=os.listdir(CleanDataDir)
    Earray=[]
    OmArray=[]
    Warray=[]
    CmeanArray=[]
    LegendList=[]
    
    SizeList=[]
    for File in FileList:
        if IsNPZ(File):
            Data = load(CleanDataDir+File)
            Elist=Data["Elist"]
            Blist=Data["Blist"]
            CmeanList=Data["Cmeanlist"]
            CmaxList=Data["Cmaxlist"]
            Om2List=Data["Om2list"]
            WList=Data["Wlist"]    
            
            [L,Nperiods,PBC]=list(Data["Parms"])
            
            Indices = Criterion(L,Nperiods,PBC,WList,Om2List)
                
            if len(Indices)>0:
                Earray.append(Elist[Indices])
                CmeanArray.append(CmeanList[Indices]/L)
                
                OmArray.append(Om2List[Indices])
                Warray.append(WList[Indices])
                

                SizeList.append([L,Nperiods,PBC])
                    
    return [Earray,CmeanArray,OmArray,Warray,SizeList]
    

    
def GroupSeries(Evec,CmeanVec,OmVec,Wvec,dW=1e-7,dOm=1e-7,sorting="omega2"):
    """Group lists after the values of the parameters omega2 and W, such that duplicate parameter sets are avoided

    args (all 1d arrays of the same length):
        Evec      : list of energy absorptions
        CmeanVec  : list of mean corr. length
        Wvec      : list of W's
        OmVec     : list of omega2's 
        
        Input should be such that [Evec[n],CmeanVec[n]] are obtained with parameters [Wvec[n],OmVec[n]] 
        
    output:
        [Wout[k],Omout[k]]: list of the Z distinct parameter sets used for the data (parameter sets are identical, if they are within [dOm,dW] from each other)
        Emean[k], Cmean[k]: mean value of E and C for the parameter set Wout[k],OmOut[k]
        Estd[k],Cstd[k] : Standard deviation of E and C  for the parameter set Wout[k],OmOut[k]
        
        Nlist : Number of data points with parameters Wout,OmOut
    """
    

    OmOut=[]
    WOut=[]
    Elist=[]
    Clist=[]

    
    for nr in range(0,len(OmVec)):
        Om=OmVec[nr]
        E=Evec[nr]
        C=CmeanVec[nr]
        W=Wvec[nr]
        
        OmMatch = abs(Om-OmOut)<dOm
        WMatch = abs(W-WOut)<dW
        
        Match=OmMatch*WMatch
        
        if sum(Match)==0:
            
            OmOut.append(Om)
            WOut.append(W)
            
            Elist.append([E])
            Clist.append([C]) 
            

        else:
            z=where(Match)[0][0]
            
            Elist[z].append(E)
            Clist[z].append(C)
    

    Elist=array(Elist)
    Clist=array(Clist)
    OmOut=array(OmOut)
    WOut=array(WOut)
    
    if sorting=="omega2":
        
        AS=argsort(OmOut)
    elif sorting=="W":
        AS=argsort(WOut)
    else:
        raise ValueError("sorting must be either 'omega2' or 'W'")
        
    OmOut=OmOut[AS]
    Elist=Elist[AS]
    Clist=Clist[AS]
    WOut=WOut[AS]
    
    Emean=array([mean(x) for x in Elist])
    Estd=array([std(x) for x in Elist])
    
    Cmean=array([mean(x) for x in Clist])
    Cstd=array([std(x) for x in Clist])    
    
    Nlist=array([len(x) for x in Elist])
    
    return [[Emean,Estd],[Cmean,Cstd],Nlist,OmOut,WOut]

def GroupData(Earray,CArray,OmArray,Warray,sorting = "omega2"):
    """ 
    Group a list of data sets in the same way as GroupSeries
    Output[n]=GroupSeries(Input[n])
    """
    Nseries=len(Earray)
    
    EmeanArray=[]
    EstdArray=[]
    CmeanArray=[]
    CstdArray=[]
    NArray=[]
    Om2Array=[]
    WArray=[]
    
    
    for ns in range(0,Nseries):
        OmVec=OmArray[ns]
        Evec=Earray[ns]
        CmeanVec=CArray[ns]
        Wvec=Warray[ns]
    
        [[Emean,Estd],[Cmean,Cstd],Nlist,Omega2,Wout]=GroupSeries(Evec,CmeanVec,OmVec,Wvec,sorting=sorting)

        EmeanArray.append(Emean)
        EstdArray.append(Estd)
        CmeanArray.append(Cmean)
        CstdArray.append(Cstd)
    
        Om2Array.append(Omega2)
        WArray.append(Wout)
        NArray.append(Nlist)
    return [[EmeanArray,EstdArray],[CmeanArray,CstdArray],NArray,Om2Array,WArray]