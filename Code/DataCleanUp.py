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
import DataAnalysisFunctions as D

OutParmList=[]
OutPathList=[]
for File in D.FileList:
    

    if D.IsNPZ(File):

        FileNameParms=D.FileNameParser(File)
        try:
            
            [L,Nperiods,PBC,omega2,W,Date,Time]=FileNameParms
        except ValueError:
        
            [L,Nperiods,PBC,Date,Time]=FileNameParms
    
                
                
        Criterion = Date > 190115 
        
        OutParm=[L,Nperiods,PBC]
        
        if Criterion:
            
            if not OutParm in OutParmList:
                OutParmList.append(OutParm)
                OutPathList.append([])
    
    
            Path = D.DataDir+File
            
            Ind=OutParmList.index(OutParm)
            
            OutPathList[Ind].append(Path)
            
    
Nparms = len(OutParmList)
def FillString(Num,width=6):
    LN=len(str(Num))
    return " "*max(0,(width-LN))+str(Num)

for nparm in range(0,Nparms):
    [L,Nperiods,PBC]=OutParmList[nparm]
    
    OutParms=[L,Nperiods,PBC]
    Elist=[]
    Blist=[]
    Cmaxlist=[]
    Cmeanlist=[]
    
    Om2list=[]
    Wlist=[]
    for Path in OutPathList[nparm]:
        try:
            
            Data=load(Path)
        except OSError:
            print(f"Bad file: {Path}. Delete?")
            x = input()
            
            if x == "y":
                os.remove(Path)
                print("    Removed file.")
            else:
                pass
#            raise OSError
        Npoints = len(list(Data["EdgeList"]))
        if Npoints ==0:
            pass #os.system(f"rm {Path}")
            
            
        else:
            E=list(Data["EdgeList"])
            B=list(Data["BulkList"])
            Cmean=list(Data["Clist"][:,0])
            Cmax=list(Data["Clist"][:,1])
            
            Parameterlist=Data["ParameterList"]
            
            [Null,Null,Null,omega2List,W]=[Parameterlist[:,z] for z in range(0,5)]
            
            omega2List=list(omega2List)
            W=list(W)
            
            Elist=Elist+E
            Blist=Blist+B
            Cmaxlist=Cmaxlist+Cmax
            Cmeanlist=Cmeanlist+Cmean
            
            Om2list=Om2list+omega2List
            Wlist=Wlist+W
 
    
    
        
    [Elist,Blist,Cmeanlist,Cmaxlist,Om2list,Wlist]=[array(x) for x in [Elist,Blist,Cmeanlist,Cmaxlist,Om2list,Wlist]]
    try:
        
        OutFileName="%d_%d_%d.npz"%(L,Nperiods,PBC)
    except TypeError:
        OutFileName="%s_%s_%s.npz"%(L,Nperiods,PBC)
    
    OutFilePath = "../DataClean/"+OutFileName
    
    savez(OutFilePath,Elist=Elist,Blist=Blist,Cmeanlist=Cmeanlist,Cmaxlist=Cmaxlist,Om2list=Om2list,Wlist=Wlist,Parms=OutParms)
    
    print(f"     L= {FillString(L,width=4)}, Nperiods = {FillString(Nperiods,width=6)}, PBC = {PBC}:  # data points: {FillString(len(Elist))}")