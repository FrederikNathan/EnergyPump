#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 11:16:53 2019

@author: frederik
"""


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
from matplotlib.pyplot import *

#DataDir = "../Data/"

L=200
Nperiods = 20000
PBC = 0

Wmax=5*2*pi
OmMax=6*2*pi
dW=0.3
dOm=0.3*1.61823


omega1=2*pi



Wrange=arange(0.0314,Wmax,dW)
OmRange=arange(0.0314,OmMax,dOm)
NW=len(Wrange)
NOm=len(OmRange)

Earray=[[[] for i in range(0,NW)] for i in range(0,NOm)]
CmaxArray=[[[] for i in range(0,NW)] for i in range(0,NOm)]
CmeanArray=[[[] for i in range(0,NW)] for i in range(0,NOm)]
Barray=[[[] for i in range(0,NW)] for i in range(0,NOm)]


DataPath=D.CleanDataDir+"%d_%d_%d.npz"%(L,Nperiods,PBC)
if not os.path.isfile(DataPath):
    raise SystemExit("No data have been generated")

Data = load(DataPath)
Elist=Data["Elist"]
Blist=Data["Blist"]
CmeanList=Data["Cmeanlist"]
CmaxList=Data["Cmaxlist"]
Om2List=Data["Om2list"]
WList=Data["Wlist"]

npoints = len(Elist)


for n in range(0,npoints):
    
    [E,B,Cmean,Cmax,omega2,W]=[X[n] for X in [Elist,Blist,CmeanList,CmaxList,Om2List,WList]]       
    Cmean=Cmean/(L/3)          
    Cmax=Cmax/(L/3)
    
    nOm=D.FindIndex(omega2,OmRange)
    nW=D.FindIndex(W,Wrange)


    if nOm!=None and nW!=None:
        Earray[nOm][nW].append(E)
        Barray[nOm][nW].append(log10(B))
        CmeanArray[nOm][nW].append(Cmean)
        CmaxArray[nOm][nW].append(Cmax)      

            


        

def FindMean(Array):
    Out=array([[nan_to_num(mean(x)+5)-5 for x in y] for y in Array])   
    Count=array([[len(x) for x in y] for y in Array])

    return [Out,Count]         


[Ogrid,Wgrid]=meshgrid(Wrange/omega1,OmRange/omega1)
vminlist=[-1.,0,0,-5,-5]
vmaxlist=[2.05,1.2,2,2]

nfig=1

TitleStr="L=%d, Nperiods = %d"%(L,Nperiods)
TitleList=[f"Edge energy absorption [$\omega_1\omega_2/2\pi$], {TitleStr}",f"Mean correlation length [$L/3$], {TitleStr}","Max correlation length [$L/3$], {TitleStr}",f"Bulk energy absorption (log10) [$\omega_1\omega_2/2\pi$], {TitleStr}"]

for Array  in [Earray,CmeanArray,CmaxArray, Barray]:
    if not ((nfig !=4  and PBC==1) or (nfig==4 and PBC==0)):
    
        figure(nfig)
        clf()
        [Mean,Count]=FindMean(Array)
        pcolormesh(Ogrid,Wgrid,Mean,cmap="gnuplot",vmin=vminlist[nfig-1],vmax=vmaxlist[nfig-1])
        colorbar()
        ylabel("$\omega_2$  $ [\omega_1]$")
        xlabel("Disorder strength  $[\omega_1]$")
        title(TitleList[nfig-1])
    
    nfig+=1
        
show()