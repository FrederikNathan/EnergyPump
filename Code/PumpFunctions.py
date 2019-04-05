#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 13:29:17 2018

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


# Core module, where the simulation is done:

# FindEnergyAbsorption is the core function. Constructs the model and simulates the model over Nperiods periods. Extracts energy absorption and correlation length.

# EnergyAbsprotionSweep is the core sweep function. Evaluates FindEnergyFunciton over a list of parameter-sets, and saves all data into a single file

def EnergyAbsorptionSweep(ParameterList,FileString):
    # Core sweep function. 
        
    # Input: ParameterList: List of parameter sets to be run over. Each element in ParameterList (parameter set) a list of the format [L,Nperiods,PBC,omega2,W].
    # FileString: File name in ../Data/ where the data are saved. 
    
    EdgeList=[]
    BulkList=[]
    Clist=[]
    DisorderList=[]

    RunID=BF.ID_gen()

    n=0
    t0=time.time()
    for Parameters in ParameterList:

        
        [E,B,C,DisorderVec]=FindEnergyAbsorption(Parameters,OutputProgress=False)
        

        Cmean=mean(C)
        Cmax=max(C)
        
        EdgeList.append(E)
        BulkList.append(B)
        Clist.append([Cmean,Cmax])
        DisorderList.append(DisorderVec)
        n+=1
      
    EdgeList=array(EdgeList)
    BulkList=array(BulkList)
    Clist=array(Clist)
        
        
    savez("../Data/%s_%s.npz"%(FileString,RunID),Clist=Clist,BulkList=BulkList,EdgeList=EdgeList,ParameterList=ParameterList)
        
        
        


def FindEnergyAbsorption(Parameters,OutputProgress=True):
    
    
    # [E,B,C,DisorderVec] = FindEnergyAbsorption([L,Nperiods,PBC,omega2,W],OutputProgress=True)
    
    # Constructs the model and simulates the model over Nperiods periods. 
    
    # L             : number of unit cells
    # Nperiods      : number of periods of mode 1 (step cycles) to run over
    # PBC           : periodic boundary conditions. If 1, periodic boundary condtitions, if 0, open boundary conditions
    # omega2        : angular frequency of ramping cycle
    # W             : Disorder strength
    
    # OutputProgress: display progress on console. 
    
    # Returns 
    
    # E:            : average energy absorption over the cycle, when initially filling the left half of the system (L leftmost sites)
    # B:            : standard deviation of energy absorption over all sites (measure of bulk absoption)
    # [Cmean,Cmax]  : mean and max correlation length of U(Nperiods *T) (i.e. the correlation length of the collumns of U)
    # DisorderVec   : Disorder potential used (can be used to reproduce results)



    [L,Nperiods,PBC,omega2,W]=Parameters
    
    L = int(L+0.1)
    Nperiods = int(Nperiods+0.1)
    PBC=int(PBC+0.1)
    
    J=2.*pi     # Tunneling strength
    omega1=2*pi              # Angular frequency of step 
    Dt= 0.25    # Duration of segment (we work in units where the driving period 
            # of mode 1 is set to 1)

    # Initialize tight binding module
    TBM.SetDimensions(OrbitalDimension=2,Dimension=1)
    
    # Initialize lattice module 
    Lattice = LM.Lattice(L,PBC=PBC) # PBC = 0 implies open boundary conditions. 
    
    # Construct Hamiltonian terms
    
    # Tunnelng term in segment 1
    T1=TBM.Hamiltonian() 
    T1[1][1,0]=-J
    T1[-1][0,1]=-J
    
    # Tunneling term in segment 3
    T2=TBM.Hamiltonian() 
    T2[1][0,1]=-J
    T2[-1][1,0]=-J
    
    # Spin operators
    sx=TBM.Hamiltonian() 
    sy=TBM.Hamiltonian()
    sz=TBM.Hamiltonian()
    
    sx[0]=TBM.SX
    sy[0]=TBM.SY
    sz[0]=TBM.SZ
    
    [SX,SY,SZ]=[sp.csc_matrix(Lattice.LatticeHamiltonian(x)) for x in [sx,sy,sz]]
    
    # Identity operator 
    Id=  Lattice.Identity()
    
    X=array([Lattice.CoordinateVec(0)]).T
    #OV = 
    DX = X.dot(ones((1,2*L)))
    DX=mod(DX-DX.T+L/2,L)-L/2
    
    # Disorder potential
    DisorderVec = W*(1-2*random((L)))
    Potential = sp.csc_matrix(Lattice.gen_OnSitePotential(DisorderVec))
    
    # Construct unitaries for first and third segment
    [H1,H2]=[Lattice.LatticeHamiltonian(x,format="csc") + Potential for x in [T1,T2]]
    
    U1 = expm(-1j*Dt*H1)    # Unitary for first segment
    U3 = expm(-1j*Dt*H2)    # Unitary for third segment
    
    X=array([Lattice.CoordinateVec(0)]).T
    DX = X.dot(ones((1,2*L)))
    DX=DX-DX.T
    
    def CouplingRange(U):
        Y = sqrt(sum(DX**2*abs(U)**2,axis=1))
  
        return Y     
    # =============================================================================
    # 3: Construct building blocks for 2nd and 4th segment unitaries
    # =============================================================================
    #
    # Unitaries for segment 2 and 4 can be computed explicitly using a rotating 
    # frame transformation. The transformation is generated by Q_{\pm }(t)  \equiv 
    # \exp{\mp i\omega_2 t \sigma_z /2}, where + corresponds to 2nd segment, and - 
    # to 4th segment. 
    # In the nth driving period, the unitary in the segment s \in {2,4} is given by 
    # 
    # U^{\pm}_n = Q_{\pm}(n T)\tilde U_{\pm} Q_{\pm}^{\dagger}(nT)
    # 
    # where \tilde U_{\pm} = Q^\dagger(sDt)e^{-i H_{\pm} Dt} Q_{\pm}((s-1)Dt)^\dagger,
    # and where H_{\pm} = JSx \pm \omega2 / 2 Sz denotes the rotating frame Hamiltonian. 
    #
    # =============================================================================
   
    # Construct Hamiltonians in rotating frame 
    Hr2 = -J*SX + omega2/2*SZ + Potential
    Hr4 = -J*SX - omega2/2*SZ + Potential
    
    # Construct unitaries 
    
    # Rotating frame unitaries over one segment
    Q2s=expm(-1j*omega2*Dt*SZ/2)
    Q4s=expm(1j*omega2*Dt*SZ/2)
    
    # Rotating frame unitaries over one period. 
    Q2p=Q2s**4# expm(-1j*omega2*T*sz)
    Q4p=Q4s**4# expm(-1j*omega2*T*sz)
    
    # Unitaries \tilde U_{\pm}
    US2 =(Q2s**2).dot(expm(-1j*Dt*Hr2)).dot(Q2s.conj())
    US4 =(Q4s**3).dot(expm(-1j*Dt*Hr4)).dot(Q4s.conj()**2)
    
    
    # =============================================================================
    # 4: Time evolve
    # =============================================================================
    if OutputProgress:
        
        print("Computing evolution operator")
    U = Id.toarray()
    t0=time.time()
    
    # 
    Q2 =  Id
    Q4 =  Id
    
    # Heisenberg picture evolution of SZ
    SZH= real(SZ.toarray().diagonal())
    SZvec= real(SZ.toarray().diagonal())
    
    # Work done on mode 1
    P = zeros((2*L,))
    
    
    def SpinZ(U):
        Out= SZvec.T.dot(abs(U)**2)
        return Out 
    
    for nT in range(0,Nperiods):
    
        for ns in range(0,4):
            t=nT+ns*Dt      
            if ns ==0 :
                dU = U1
                
            if ns ==1:
                dU = Q2.dot(US2.dot(Q2.conj()))
    
            if ns ==2:
                dU = U3
    
            if ns ==3:
                dU = Q4.dot(US4.dot(Q4.conj()))            
    
            # Time-evolve
            U = dU.dot(U)
            
            # Evolve rotating frame unitaries 
            Q2 = Q2.dot(Q2p)
            Q4 = Q4.dot(Q4p)
            
            # Compute evolution of SZ[i,i]
    
    
            if ns ==1 or ns ==3 :
                        
                SZ0=SZH
                SZH=SpinZ(U)
                if ns==1:
                    s=1
                else:
                    s=-1
                    
                dP=s*(SZH-SZ0)*0.5
                
                P=P+dP 
                           
    
        if nT%(max(1,int(Nperiods/20)))==0 and OutputProgress:
            print(f"   at period {nT}/{Nperiods}. Time spent: {time.time()-t0:.2g} s ")
            
    P=P/Nperiods 
    
    EdgeEnergyAbsorption = sum(P[:L])
    BulkEnergyAbsorption = std(P)
    
    C = sort(CouplingRange(U))

    return [EdgeEnergyAbsorption,BulkEnergyAbsorption,C,DisorderVec]




    
    
