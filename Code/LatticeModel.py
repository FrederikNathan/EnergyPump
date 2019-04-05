#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 11:10:57 2018

@author: frederik
"""

import sys
import os 
import os.path as op

from scipy import *
from scipy.linalg import *
from numpy import ufunc as ufunc
import scipy.sparse as sp
from numpy.random import *
from scipy.misc import factorial as factorial

import datetime
import logging as l
import time

import TightBindingModel as TBM 

#import scipy.ufunc as ufunc

Module=sys.modules[__name__]

SX = array(([[0,1],[1,0]]),dtype=complex)
SY = array(([[0,-1j],[1j,0]]),dtype=complex)
SZ = array(([[1,0],[0,-1]]),dtype=complex)
I2 = array(([[1,0],[0,1]]),dtype=complex)

def CheckDimensionSet():
    if not TBM.DimensionSet:
        raise ValueError("Dimension and orbital dimension must be set using SetDimension")
        
class Lattice:
    """Lattice object -- contains all lattice-specific functions and objects
    Can be combined with Clean Hamiltonian to generate a Hamiltonian matrix"""
    def __init__(self,LatticeDimension,PBC=False):
        
        if type(LatticeDimension)==int:
            LatticeDimension=(LatticeDimension,)
            
        CheckDimensionSet()
            
        self.CheckLatticeDimension(LatticeDimension)
        self.__PBC=PBC
        self.__LatticeDimension=LatticeDimension
        self.__Hdimension=TBM.OrbitalDimension*prod(LatticeDimension)
        
            
        self.__UnitCells=prod(LatticeDimension)
        
        self.__Ilist=self.gen_Ilist()
        self.__Tlist=self.gen_Tlist()
#        
#        global Hdimension 
#        
#        Module.Hdimension = self.Hdimension
#        Module.LatticeDimension = LatticeDimension
#        Module.OrbitalDimension = OrbitalDimension
#        Module.PBC = PBC

    def Ilist(self):
        return self.__Ilist
    
    def Tlist(self):
        return self.__Tlist

    def CheckLatticeDimension(self,Ldim):
        if not len(Ldim)==TBM.Dimension:
            raise ValueError("Lattice dimension must be given by n integers, where n=%d is the model's dimensionality "%TBM.Dimension)
        
        
    def get_index(self,Coords,Orbitals):
        if amax(Orbitals) >= self.OrbitalDimension:
            raise ValueError("Orbital index exceeds orbital dimension")
            
        Index = Orbitals
        Q=self.OrbitalDimension
        
        for d in range(0,TBM.Dimension):
            Index = Index + Coords[d]*Q    
            Q=Q*self.LatticeSize[d]
        
        return Index 
    
    def gen_Tlist(self   ):
        """ Get translation operators on lattice (should be directly multiplied together with identity matrices to become operators on full hilbert space of lattice)"""
        Tlist=[]
        for d in range(0,TBM.Dimension):
            D=self.__LatticeDimension[d]
            
            
            Data=ones(D)

            if self.__PBC==False:
                Data[D-1]=0
            elif self.__PBC==True:
                pass
            elif type(self.__PBC)==tuple:
                Data[D-1]=self.__PBC[d]
                
                
            Cols=arange(D)
            Rows=(Cols+1)%D
    # 
            T=sp.csr_matrix(sp.coo_matrix((Data,(Rows,Cols)),dtype=complex))
            
            Tlist.append(T)
        
        return Tlist

    def gen_Ilist(self):
        """ Get identity operators on lattice (should be directly multiplied together to become full identity operator)"""

        Ilist=[] #csr_matrix(eye(self.OrbitalDimension),dtype=complex)]
        
        for d in range(0,TBM.Dimension):
            D=self.__LatticeDimension[d]
            Ilist.append(sp.csr_matrix(eye(D),dtype=complex))
        return Ilist
    
    def Identity(self,format="csr"):
        I = sp.eye(self.__Hdimension,format=format)
        return I
    
    def CoordinateVec(self,Dim):
        """ Get list of Dim-coordinates corresponding to indices""" 
        if Dim + 1 > TBM.Dimension:
            raise IndexError("Direction argument must match the lattice dimension")
            
        LatticeIndices = arange(0,self.__Hdimension)//TBM.OrbitalDimension
        Vec=(LatticeIndices//prod(self.__LatticeDimension[:Dim]))%self.__LatticeDimension[Dim]
        
        return Vec
    
    def CoordinateOperator(self,Dim,format="csr"):
        """ Get operator corresponding to the Dim-Position operator"""
        CoordinateVec=self.CoordinateVec(Dim)
        
        Mat = sp.eye(self.__Hdimension,format=format)
        Mat.setdiag(CoordinateVec)
        return Mat
    
        

    
    def gen_OnSitePotential(self,Potential):
        
        """ 
        Generates a Hamiltonian corresponding to on-site potential
        
        The on-site potential A must be an ndarray, where A[nx,ny,..] indicates the on-site potential on site [nx,ny,...]
        """
        
        # Checking for right input"
        
        if not type (Potential)==ndarray or shape(Potential)!=tuple(self.__LatticeDimension):
            print(f"Shape of potential : {shape(Potential)}")
            print(f"Lattice dimension  : {self.__LatticeDimension}")
            raise TypeError("Shape of potential must be ndarray and match physical Hamiltonian")
            
        Vec = ravel(Potential,order="F")
        
        Data=Vec
        Rows=arange(0,self.__UnitCells)
        Cols=arange(0,self.__UnitCells)
        LatticeMatrix = sp.csr_matrix((Data,(Rows,Cols)),dtype=complex)
    
        I = sp.csr_matrix(eye(TBM.OrbitalDimension),dtype=complex)
    
        return sp.kron(LatticeMatrix,I,format="csr")    
    

    def get_density(self,Vec):
        
        Rho=abs(Vec)**2
        
        if len(shape(Rho))>1:
            
            Shape=tuple([TBM.OrbitalDimension]+list(self.__LatticeDimension)+[shape(Rho)[1]])
            
        else:
            Shape=tuple([TBM.OrbitalDimension]+list(self.__LatticeDimension))   
        RhoMat=reshape(Rho,Shape,order="F")
        
        RhoLattice = sum(RhoMat,axis=0)
        return RhoLattice

    
    
    def VecToMatrix(self,Vec):
        Shape=tuple([TBM.OrbitalDimension]+list(self.__LatticeDimension))
        
        return reshape(Vec,Shape,order="F")

    def get_Identity(self,format="csr"):
        
       return sp.eye(self.__Hdimension,format=format,dtype=complex)

    

        
    def MatrixToVec(Mat):
        return ravel(Mat,order="F")    
    

    
    def get_coords(index):
        
        orbital = index%TBM.OrbitalDimension
        latticeindex=index//TBM.OrbitalDimension
        
        coordlist=[]
        for d in range(0,Dimension):
            coord=(latticeindex//prod(self.__LatticeDimension[:d]))%self.__LatticeDimension[d]
            coord=cottord.astype(int)
            
            coordlist.append(coord)
            
        return [coordlist,orbital]


    def LatticeHamiltonian(self,Hamiltonian,format="csr"):
        """ Generates real space Hamiltonian matrix corresponding to lattice. 
        orbital d, position r=[x,y,z] corresponds to index d+x D_orbital + y D_orbital Lx + z D_orbital Lx Ly"""
#    #    CheckLattice()
#        print(type(Hamiltonian))
#        if not type(Hamiltonian)==TBM.Hamiltonian:
#            print(type(Hamiltonian))
#            raise ValueError("Hamiltonian must be BZO Hamiltonian")
        
        Hdimension = self.__Hdimension
        Dimension=TBM.Dimension   
        Tlist=self.Tlist()
        
        NNZ=Hamiltonian.NNZ()
        IndList=Hamiltonian.IndList()
        
        Matrix = sp.csr_matrix((Hdimension,Hdimension),dtype=complex).asformat(format)
        for Ind in IndList:
            OrbitalMat = Hamiltonian[Ind]
    
            Out = sp.csr_matrix(OrbitalMat,dtype=complex)
    
            for d in range(0,TBM.Dimension):
                if Ind[d]>=0:
                    power = Ind[d]
                    Mat = Tlist[d]
                if Ind[d]<0:
                    power = - Ind[d]
                    Mat = Tlist[d].conj().T
                    
                t = pow(Mat,power)
                
    
                
                Out = sp.kron(t,Out,format=format)
    
            
    
    
            Matrix = Matrix + Out
        
        return Matrix



