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


Module=sys.modules[__name__]

SX = array(([[0,1],[1,0]]),dtype=complex)
SY = array(([[0,-1j],[1j,0]]),dtype=complex)
SZ = array(([[1,0],[0,-1]]),dtype=complex)
I2 = array(([[1,0],[0,1]]),dtype=complex)
#
#Dimension = 3
#OribtalDimension= 2 

DimensionSet=False

def SetDimensions(OrbitalDimension=None,Dimension=None):
    Module.Dimension=Dimension
    Module.OrbitalDimension=OrbitalDimension
    Module.DimensionSet=True
    


def CheckDimensionSet():
    if not DimensionSet:
        raise ValueError("Dimension and orbital dimension must be set using SetDimension")
        
        
def CheckLattice():
    if not "Lattice" in globals():
        raise NameError("Lattice object must be specified")
    elif not type(Lattice)==LatticeObject:
        print(type(Lattice))
        raise NameError("Lattice must be a lattice object.")
        
    
        
class BZO:
    """Base BZ field object. 
    
    Since a field on the BZ is a periodic function of 
    crystal momentum, it can be decomposed into discrete Fourier Harmonics. 
    Hamiltonians have particularly simple harmonics. 
    
    The BZO field saves the object through its Fourier components (saved in ObjList),
    with the corresponding "wavevector" saved in IndList. 
   
    All standard operations, as well as scalar addition are defined for BZO's.
    This included multiplication and addition of two BZOs (resulting in a third 
    BZO)
    
    When BZO is called at a crystal momentum, or an array or list of crystal momenta,
    an array is returned with the values of the field at the k-points listed in the array. 
    10 x  efficiency can be achieved when calling the BZO function with the vector span-format (see __call__).
    
    When list element is referred to, returns corresponding harmonic, such that BZO[n,m,k] returns the 
    (n,m,k)th Fourier Harmonic """
       
    
    """Also contains numlist, which sorts the indices, and is used for rapidly finding elements in 
    indlist and objlist. In this list, indices ind IndList are translated to integers, using
    the __IndexCounter variable. __IndexCounter is not fixed, but is by default set to 1e6.
    Elements in IndList and ObjList are sorted according to their values in NumList. """
    
    
    """ If list element (a,b) is empty, F[a,b] returns BZO.__ZeroObject. This object acts as the 
    zero matrix in all ways, except that a new object and index is added to the IndList and ObjList, if 
    F[a,b] is set to a nonzero value. 
    I.e. 
    
    print(4*F[a,b]) returns  [[0,0]
                            [0,0]]
    
    But if we set F[a,b] = X
    
    (a,b) is added to indlist, with corresponding ObjList element giben by  x
    such that 
    
    F[a,b]=x
    
    """ 
    
    def __init__(self,*args,shape=None,dtype=complex):# shape,dtype=complex):

        CheckDimensionSet()

        ArgV=[arg for arg in args]
        
        if len(ArgV)==2:
            if not shape==None:
                print("Warning: shape is not considered, but follows from field argument")
                
            [IndList,ObjList]=ArgV

            self.__Dimension = Dimension
            
            self.__shape=array(ObjList[0]).shape
            self.__dtype=ObjList[0].dtype
            
            self.__Zero = self.__ZeroObject(self)               
            self.SetLists(IndList,ObjList)
            
        elif len(ArgV)==1:
            if not shape==None:
                print("Warning: shape is not considered, but follows from field argument")
                
            Obj = generate_field(ArgV[0])

            self.__shape=Obj.shape()
            self.__dtype=Obj.dtype()
            
            self.__Dimension = Dimension
            self.SetLists(Obj.IndList(),Obj.ObjList())

            self.__Zero = Obj.__Zero


            
            
            
        elif len(ArgV)==0:
            if not (type(shape)==tuple or type(shape)==int):
                raise ValueError("Shape must be specified, and given by an integer or a tuple")

            self.__Dimension=Dimension
            
            self.SetLists([],[])
                        
            self.__shape=shape
            self.__dtype=dtype
            
            self.__Zero = self.__ZeroObject(self)
            
        else:
            raise ValueError("BZO can only take a single array of other BZO's as argument")
  
    
    def SetLists(self,IndList,ObjList):
        """ Set IndList and ObjList for object"""
        
        if not len(IndList)==len(ObjList):
            raise ValueError("IndList and ObjList must be of the same length")
            
        if not prod([type(x)==tuple for x in IndList]):
            IndFormatError=1
        elif not prod([len(x)==3 for x in IndList]):
            IndFormatError=1
        
        IndFormatError=0
        
        if IndFormatError:
            raise ValueError("IndList must be list of 3-tuples")
            
        if not prod([shape(x) == self.__shape for x in ObjList])==1:
            print([shape(x) for x in ObjList])
            raise ValueError("ObjList must contain arrays of the same shape as the BZO (shape %s)"%str(self.__shape))
            
            
        self.__IndList=IndList
        self.__ObjList=ObjList
        
        self.__Set_NumList()    
        self.__SortLists()
        
        
    def __Set_NumList(self,Num=None):
        """ Compute NumList for object"""
        DefaultIndexCounter=1e5
        
        if Num == None:
            if len(self.__IndList)==0:
                IndexCounter = DefaultIndexCounter
                
            else:
                IndexMax = max([max([abs(x) for x in Ind]) for Ind in self.__IndList])
                IndexCounter = int(max(DefaultIndexCounter,IndexMax)+0.1)

        else:
            IndexCounter = Num
            
        self.__IndexCounter=IndexCounter
        
        self.__NumList=[self.__IndToNum(Ind) for Ind in self.__IndList]
        

    def __IndToNum(self,*Ind):
        """ Recipe for converting Index in IndList to integer. Used for sorting IndList and ObjList (see NumList),
        and for making it __FindIndex go faster. Uses the variable __IndexCounter, which by default is set to 1e6, but changes if larger values are required """
        
        if self.__Dimension>1:
            return sum((Ind[0][n]*self.__IndexCounter**(self.__Dimension-n-1) for n in range(0,self.__Dimension)))     
        
        else:
            try:
                return Ind[0][0]*self.__IndexCounter**(self.__Dimension-1)     
            except TypeError:
                return Ind[0]*self.__IndexCounter**(self.__Dimension-1)




    def __CheckIndices(self,Index):
        """ Check that index has the right format"""
        if not len(Index) == self.__Dimension:

            raise IndexError("Field is indiced by %d integers. Index given is %s"%(self.__Dimension,str(Index)))
        
    def __FindIndex(self,*Index):
        """ Find the list index corresponding to *Index. If *Index is not in 
        IndList, return -1-n where n marks the list slot to the right of where *Index 
        would go if it were to be inserted"""
        
 
        self.__CheckIndices(*Index)
        
        Num=self.__IndToNum(*Index)
        listindex= searchsorted(self.__NumList,Num)



        if listindex<self.NNZ():

            if self.__NumList[listindex]==Num:
                
                return listindex
            
            else:
                
                return -1-listindex
            
        else:
            return -1-listindex
            

    def __getitem__(self,*Index):

        if type(*Index)==int:
            Index = (Index,)
            

        n = self.__FindIndex(*Index)
        
        if n>=0 :

            return self.__ObjList[n]
        
        else:
            

            self.__Zero.CalledIndex=tuple(*Index)  

            return self.__Zero
        

    class __ZeroObject:
        def __init__(self,bzo):

             self.__shape=bzo.shape()
             self.__Mat = zeros(self.__shape,dtype=complex)
             self.__BZO = bzo
             self.CalledIndex=None
             
             
        def __getitem__(self,Index):
             return self.__Mat[Index]*1
         
        def __setitem__(self,Index,Value):
            NewMat = 1*self.__Mat

            NewMat[Index]=Value
            
            self.__BZO[self.CalledIndex]=NewMat
            
    
        def __str__(self):
            return str(self.__Mat)


        def __repr__(self):
            return self.__str__()
             
        
        def __add__(self,x):
            return self.__Mat + x
        
        def __radd__(self,x):
            return self+x
        
        def __mul__(self,x):
            return self+0
        
        def __rmul__(self,x):
            return self * x 
        
        
        def __sub__(self,x):
            return self+(-x)
        
        def __rsub__(self,x):
            return (-1)*self + x 
        
        def truediv(self,x):
            return (1/x)*self
        

        
    def __setitem__(self,Index,Value):
        # Check that the value has the right format 

        if type(Index)==int:
            Index = (Index,)
            
            
        Value=array(Value)
        if (not type(Value)==ndarray) or (not shape(Value)==self.__shape):
            raise ValueError("Assigned value must be ndarray of shape %s"%str(self.__shape))
                  
        Q = 3*max([abs(n) for n in Index])
        
        if Q>self.__IndexCounter:
            self.__Set_NumList(Num=Q)
            
        n = self.__FindIndex(Index)
        
        
        if n >=0:
            self.__ObjList[n]=Value
            
        else:
            newIndex = -n-1
            
            self.__IndList.insert(newIndex,Index)
            self.__ObjList.insert(newIndex,Value)
            self.__NumList.insert(newIndex,self.__IndToNum(Index))
            

    


        
    def __SortLists(self):
        """ Sorts indlist, objlist and numlist according to their NumList values
        Elements are first sorted according to the first index of Ind, then the 
        second index and so on. """ 

        
        AS=argsort(self.__NumList)

        self.__IndList=[self.__IndList[i] for i in AS]#list(self.__IndList[AS])
        self.__ObjList=[self.__ObjList[i] for i in AS]#list(self.__IndList[AS])
        self.__NumList=[self.__NumList[i] for i in AS]
        
    def NNZ(self):
        """ Returns number of nonzero elements in ObjList"""
        return len(self.__IndList)
    
    def ObjList(self):
        """ Returns ObjList"""
        return self.__ObjList[0:]*1
    
    def IndList(self):
        """ Returns IndList"""
        
        return self.__IndList[0:]*1       #Use 0: to ensure that returned object is plain data, not a reference to self.__IndList
    
    def __mul__(self,y):     
        """ Multiplication. If F1(k) is multiplied with scalar lambda, returns lambda * F(k)
        If mulitplied with BZO of same shape, returns Out(k) = F1(k)*F2(k)""" 

        # BZO mulitplication
        if type(y)==type(self):
            Out = self._CreateSameType()
            
            for Ind1 in self.IndList():
                Obj1=self[Ind1]
                for Ind2 in y.IndList():
                    Obj2=y[Ind2]
                    
                    Ind3 = tuple(add(Ind1,Ind2))
    
                    Out[Ind3] += Obj1*Obj2
     
        # Scalar multiplicatin
        else:

            Out = self._CreateSameType()

            Out.SetLists(self.IndList(),[y*x for x in self.__ObjList])

            # Multiplication with item of its own type
            
                
            
                    
        
        return Out
    
        
    def __rmul__(self,x):
        return self*x

    
    def __truediv__(self,y):
        return self*(1/y)
    

    def __add__(self,Obj):

        if not type(self)==type(Obj):
            raise ValueError("Two added objects must be of the same type. Type of arguments are "+str(self.__class__.__name__)+", "+str(Obj.__class__.__name__))

        Out = self._CreateSameType()
    
    
        ObjListOut = self.ObjList()
        IndListOut = self.IndList()
        

        IndList2 = Obj.IndList()
        ObjList2 = Obj.ObjList()       


        if self.__IndexCounter == Obj.__IndexCounter:
            N1=1*self.__NumList
            N2=1*Obj.__NumList
            
            NNZ = len(N1)
            NNZ2 = len(N2)



            z1 = 0
            z2=0
           
            while z2 < NNZ2 and z1 < NNZ:
                Num2 = N2[z2]

                
                if Num2 ==N1[z1]:
                    ObjListOut[1*z1]=1*(ObjListOut[z1]+ObjList2[z2])
                    z2+=1


                elif Num2 < N1[z1]:
                    ObjListOut.insert(z1,ObjList2[z2])
                    IndListOut.insert(z1,IndList2[z2])
                    N1.insert(z1,Num2)
                    
                    z2+=1
                    NNZ+=1
                    z1+=1

                elif Num2 > N1[z1]:
                        
                    z1+=1
                        
            ObjListOut = ObjListOut + ObjList2[z2:]
            IndListOut = IndListOut + IndList2[z2:]

                
            Out.SetLists(IndListOut,ObjListOut)
            Out.CleanUp()               
            return Out
            

        else:

            Out.SetLists(IndListOut,ObjListOut)
            
            for Ind in IndList2:

                Out[Ind]+=Obj[Ind]

            Out.CleanUp()
            return Out
    
    def __radd__(self,Obj):

        return self + Obj
    
    def _CreateSameType(self):
        
        return BZO(shape=self.__shape,dtype=self.__dtype)

    def __sub__(self,Obj):
        
        return (-1)*Obj + self

    def __rsub__(self,Obj):
        
        return self-Obj
  
    def __InputInterpreter(self,*args):
        """Interprets input from function call. Input can be of 3 types: VectorSpan, List, or Array"""
        
        if len(args)==self.__Dimension:
            VectorLists = [array(x,ndmin=1) for x in args]
        
            Format="VectorSpan"
            
            return [VectorLists,Format]
        
        elif len(args)==1:
            
            K = args[0]
            
            if type(K)==list:
                Karray=array(args[0]).T
                
                Format="List"
            else:
                Karray=K
                
                Format="Array"
                
            
            if (not type(Karray)==ndarray) or (not shape(Karray)[0]==self.__Dimension):
                raise IndexError("Argument must be ndarray of dimension %d x *"%self.__Dimension )
                
            return [Karray,Format]
        else:
            raise ValueError("Did not understand input")
            
    
    def __FindIndexSpan(self):
        
        self.__IndexSpanList=[]
        self.__IndexSpanPointerList=[]
        
        for d in range(0,self.__Dimension):
            
            List = list(set([Ind[d] for Ind in self.__IndList]))
            List.sort()
            
            self.__IndexSpanList.append(List)
            
            IndexSpanPointer  = [searchsorted(List,Ind[d]) for Ind in self.__IndList]
            
            self.__IndexSpanPointerList.append(IndexSpanPointer)
        
        
        
    def __call__(self,*args):
                     
        """ Evaluate function at kvec (kvec is a momentum or array of momenta). 
        Returns value, or array of values"""
       
        [Input,Format]=self.__InputInterpreter(*args)
      
        
        if Format=="List" or Format=="Array":
            Karray =Input
            
                
            OutShape = shape(Karray)[1:]
            
            Karray_dimension=len(OutShape)
            
            OutMatShape = self.__shape+OutShape
            OutMat = zeros(OutMatShape,dtype=self.__dtype)




            ## List argument to be passed to Einstein summation
            EinsteinList = [x for x in range(0,Karray_dimension+1)]
            
            ## Generate lists with exponential vectors
   
    
            for n in range(0,self.NNZ()):
    
                Ind = array(self._BZO__IndList[n])
                C   = self._BZO__ObjList[n]           
                
                   
                PhaseMat = exp(-1j*einsum(Ind,[0],Karray,EinsteinList))
                
                if not self.__shape ==():
                    
                    dO = multiply.outer(C,PhaseMat)
                else :
                    dO = C*PhaseMat
                    
                OutMat = OutMat + dO
                
    
            
    
            if Format=="Array":
                return OutMat        
            
            elif Format=="List":
                return list(OutMat)
            
        elif Format=="VectorSpan":
            
            Vectors= Input
            
            OutShape = tuple([len(v) for v in Vectors])
            
            
            OutDimension=len(OutShape)
            
            OutMatShape = self.__shape+OutShape
            OutMat = zeros(OutMatShape,dtype=self.__dtype)
            
            
            self.__FindIndexSpan()
            
            VectorList=[]
            MultList=[]
            
            Module.OutArrayList= []
            for d in range(0,self.__Dimension):
                IS = self.__IndexSpanList[d]
                V = Vectors[d]
                MultList.append( [exp(-1j*V*Index) for Index in IS])

                
            Module.OutArrayList = [zeros(OutShape[-d:],dtype=complex) for d in range(2,self.__Dimension+1)]
            dO = zeros(OutMatShape,dtype=complex)
              

            def Multiply(Args):
                Nargs = len(Args)
                if Nargs>1:

                    return multiply.outer(Args[0],Multiply(Args[1:]),out=Module.OutArrayList[Nargs-2])
                else:
                    return Args[0]
#                
                
            
            for n in range(0,self.NNZ()):

                C   = self._BZO__ObjList[n]           
                
                
                VecList=[MultList[d][self.__IndexSpanPointerList[d][n]] for d in range(0,self.__Dimension)]

                    
                PhaseMat = Multiply(VecList)

                if not self.__shape ==():
                    
                    multiply.outer(C,PhaseMat,out = dO)
                    
                else :
                    multiply(C,PhaseMat,out = dO)
                    
                OutMat = OutMat + dO
                
#                print(PhaseMat)
                                
            return OutMat
        
    
    def Gradient(self):
        """ Returns gradient \partial _{k_i} F(k) . 
        BZO.Gradient() Returns BZO field of dimension (Dimension,)+shape(BZO), with 
        BZO.Gradient()(k)[i,:,..] =\partial _{k_i} F(k)"""
        
        OutList  = []
        for dim in range(0,Dimension):
            X = self._CreateSameType()
            
            for Ind in self.IndList():
                
                Factor = -1j * Ind[dim]
                
                X[Ind] = Factor * self[Ind]
                
            OutList.append(X)
            
        return BZO(array(OutList))

    
    def shape(self):
        """ Returns shape of BZO"""
        return self.__shape
    
    
    def dtype(self):
        """Returns data type of BZO"""
        return self.__dtype 
    
    def __repr__(self):
        
        if self.NNZ()>0 and self.NNZ ()< 20:
            
            Str = "Type: %s"%str(self.__class__.__name__)+"\n\n"
            
            for Ind in self.__IndList:
                Str += str(Ind)+": \n \n"
                Str += self[Ind].__str__()
                Str += "\n \n \n"
            
            
        elif self.NNZ()==0:
            Str = "%s object with all zeros"%str(self.__class__.__name__)
            
        elif self.NNZ()>20:
            Str = "%s of shape %s with %d Matrices (too long to show here) " %(self.__class__.__name__,str(self.shape()),self.NNZ())

        return Str 
    
    def __str__(self):
        return self.__repr__()
    
    
    def conj(self):
        """ Returns conjugate transpose of BZO"""
        
        Out = self._CreateSameType()
        
        for Ind in self.IndList():
            OutInd = tuple(-x for x in Ind)
            
            Out[OutInd]=self[Ind].conj().T
            
        return Out 
    
    def slice(self,*Indices):
        """ Return BZO such that, with Y=BZO.slice(Indices),  Y(k) = BZO(k)(Indices)."""
        
        Ind = tuple(Indices)


        try:
            
            OutShape=shape((1*self[(0,)*Dimension])[Indices])
        except:
            raise IndexError("Wrong format for indices")
            
        Out = BZO(shape=OutShape)
        
        for Ind in self.IndList():

            Out[Ind]=array(self[Ind][Indices])
            
        Out.CleanUp()
        
        return Out

    def __delitem__(self,Ind):
    
        n=self.__FindIndex(Ind)
        
        if not n==None :

            del self.__IndList[n]
            del self.__ObjList[n]
            del self.__NumList[n]
        
    def CleanUp(self):
        """ Remove negligible elements in IndList and ObjList. By default, 
        elements are discared, if they have max-norm less than 1e-10
        """
        for Ind in self.IndList():
            if amax(abs(self[Ind]))<1e-10:
                del self[Ind]                
    

    
    def Norm(self):
        """ Returns norm \int d^D k \Tr (F^\dagger(k) F(k))"""
        
        return sqrt(sum([sum(abs(x)**2) for x in self.__ObjList]))
       
def Karray(*args):
    """TBA      """
    return None 

def generate_field(FieldArray):
    """Combine field from array of fields"""
    shape0=shape(FieldArray)
    
    FlatArray = FieldArray.flatten()
    N_elements = len(FlatArray)

    try:
            
        shape1=FlatArray[0].shape()
        type1=type(FlatArray[0])
        
        
        SameType = prod([type(x)==type1 for x in FlatArray])
        SameShape = prod([x.shape()==shape1 for x in FlatArray])
        
        
        if not (SameType and SameShape):
            raise ValueError("Argument must be array of objects of the same type and shape")
    except:
        raise ValueError("Argument must be array of objects of the same type and shape")
        
    
    OutShape = tuple(shape1)+tuple(shape0)
    
    Out = BZO(shape=OutShape)
        

    
    
    ### Construct IndList 
    

    SliceTuple=tuple(slice(x) for x in shape1)
    
    for n in range(0,N_elements):
        field = FlatArray[n]
        
        IndList=field.IndList()
              
        
        for Ind in IndList:
                
            MatrixIndex = unravel_index(n,shape0)
            
            Out[Ind][SliceTuple+MatrixIndex]+=field[Ind]

    return Out

       
class Hamiltonian(BZO):
 
    """Physical Hamiltonian object. Constains information about the translationally invariant part of the Hamiltonian"""
    
    def __init__(self,*args):
        BZO.__init__(self,shape=(OrbitalDimension,OrbitalDimension))        

        
        ArgV=[arg for arg in args]  
        
        if len(ArgV)==1:
  
            Obj = ArgV[0]


            if (not type(Obj)==Module.BZO) or (not Obj.shape() ==(OrbitalDimension,OrbitalDimension)):
                raise TypeError("Argument must be a field of dimension %d x %d "%(OrbitalDimension,OrbitalDimension))
            
            self.SetLists(Obj.IndList(),Obj.ObjList())
            
            if not self.IsHermitian():
                print("Warning: Hamiltonian is non-hermitian")
            
        elif len(ArgV)==0:
      
            self.SetLists([],[])
                        
        else:
            raise ValueError("BZO can only take a single array of other BZO's as argument")


        

    def Bands(self,Karray):
        """ Compute energy band structure at k-points specified by Karray""" 
        
        B = self(Karray)
        
        OutShape = (OrbitalDimension,)+shape(Karray)[1:]
        Nk = int(prod(shape(Karray)[1:])+0.1)
        
        BFlatShape=(OrbitalDimension,OrbitalDimension,Nk)

        Bflat = reshape(B,BFlatShape)

        Bands = array([eigh(Bflat[:,:,n].T,eigvals_only=True) for n in range(0,Nk)]).T
        Bands = reshape(Bands,OutShape)
        
        return Bands
    
    
    def Trace(self):
        """Compute trace of self. Returns Scalar object"""
        
        Out = Scalar()
        
        for Ind in self.IndList():
            Out[Ind]=trace(self[Ind])
            
        
        return Out 
    
    def _CreateSameType(self):
        return Hamiltonian()
    
    def IsHermitian(self):
        """Determine whether Hamiltonian is Hermitian"""
        
        Hermitian=True
        for Ind in self.IndList():
            Q=tuple(-x for x in Ind)
            
            X = self[Ind].conj().T-self[Q]

            A=amax(list(abs(X.flatten())))

            
            if A > 1e-9:
                Hermitian=False
                
        return Hermitian
    

    
class Scalar(BZO):
    """ Sxcalar-valued function on the BZ"""
    
    def __init__(self,*args):  
        if len(args)==2:
            BZO.__init__(self,*args)
        elif len(args)==0:
            BZO.__init__(self,shape=())


    def __getitem__(self,Index):
        return BZO.__getitem__(self,Index)*1

    def __setitem__(self,Index,Value):
        BZO.__setitem__(self,Index,array(Value))
        
        
    def __add__(self,y):

        if not (type(y)==Module.Scalar or type(y)==Module.BZO or type(y)==Module.Hamiltonian):
                
            Out = Scalar(self.IndList(),self.ObjList())
#            
#            Out._BZO__ObjList=self.ObjList()
#            Out._BZO__IndList=self.IndList()
#            

            Out[(0,)*Dimension]=Out[(0,)*Dimension]*1+y

        else:
            Out = BZO.__add__(self,y)
            
        return Out   
        
    
    def __radd__(self,y):
        return self+y
    
    def __sub__(self,y):

        return self + (-1*y)
    
    def __rsub__(self,y):
        return self-y
    
    def __mul__(self,y):

        if type(y) == ndarray:
                
            Out = BZO(shape=shape(y))
            
            for Ind in self._BZO__IndList:
                
                Out[Ind] = self[Ind]*y
                
            
            return Out

        else:

            return BZO.__mul__(self,y)
                
            
        
        
    def __rmul__(self,y):
        return self*y
    
    def _CreateSameType(self):
        return Scalar()


def DataToScalar(Data,Cutoff=1e-3):
    """ Turn data into scalar object. Discard scalar elements smaller than cutoff"""     
    
    import numpy.fft as fft
    
    # Renormalizing cutoff
    Cutoff = Cutoff * sqrt(sum(abs(Data)**2)/prod(shape(Data)))
    
    DataDimension = shape(Data)
    FourierData=fft.fftn(Data)/prod(DataDimension)
    
    RawIndList=where(abs(FourierData)>Cutoff)
    
    Nind = len(RawIndList[0])
    
    def FindIndex(Index):
        Out =[]
        for n in range(0,len(Index)):# in Index:
            x = Index[n]
            Dim=DataDimension[n]
            y = int((x+Dim/2)%Dim - Dim/2 )
    
            Out.append(-y)
        
        return tuple(Out)
    
    IndList= []    
    ObjList = []
    
    print("Converting data to Scalar object. Number of nonzero elements: %d"%Nind)
    for n in range(0,Nind):

        Raw_index = tuple(int(A[n]+0.1) for A in RawIndList)
        
        Ind = FindIndex(Raw_index)
        IndList.append(Ind)
        
        ObjList.append(FourierData[Raw_index])
    
    
    
    Out=Scalar(IndList,ObjList)
    print("")
    return Out
        
def GetTranslationOperators():
    OutList=[]
    for d in range(0,Dimension):
        X = Scalar()
        
        Dir = tuple(d*[0]+[1]+[0]*(Dimension-d-1))
        
        X[Dir]=1 
        
        OutList.append(X)
        
    return OutList


    
        
class LatticeObject:
    """Lattice object -- contains all lattice-specific functions and objects
    Can be combined with Clean Hamiltonian to generate a Hamiltonian matrix"""
    def __init__(self,LatticeDimension,PBC=False):
        
        CheckEnvironmentVars()
        
        if type(LatticeDimension)==int:
            LatticeDimension=tuple([LatticeDimension])
        
        # Checking for right input
        if not shape(LatticeDimension)[0]==Dimension:
            raise ValueError("Dimension of lattice sites must match dimension of syste")
                    
        self.OrbitalDimension=OrbitalDimension
        self.PBC=PBC
        self.Dimension=Dimension
        self.Hdimension=self.OrbitalDimension*prod(LatticeDimension)
        self.LatticeDimension=LatticeDimension
        self.UnitCells=prod(LatticeDimension)
        
        self.Ilist=self.gen_Ilist()
        self.Tlist=self.gen_Tlist()
#        
#        global Hdimension 
#        
        Module.Hdimension = self.Hdimension
        Module.LatticeDimension = LatticeDimension
        Module.OrbitalDimension = OrbitalDimension
        Module.PBC = PBC
    
    def get_index(self,Coords,Orbitals):
        if amax(Orbitals) >= self.OrbitalDimension:
            raise ValueError("Orbital index exceeds orbital dimension")
            
        Index = Orbitals
        Q=self.OrbitalDimension
        
        for d in range(0,self.Dimension):
            Index = Index + Coords[d]*Q    
            Q=Q*self.LatticeSize[d]
        
        return Index 
    
    def gen_Tlist(self   ):
        Tlist=[]
        for d in range(0,self.Dimension):
            D=self.LatticeDimension[d]
            
            
            Data=ones(D)

            if not self.PBC:
                Data[D-1]=0
                
                
            Cols=arange(D)
            Rows=(Cols+1)%D
    # 
            T=csr_matrix(coo_matrix((Data,(Rows,Cols)),dtype=complex))
            
            Tlist.append(T)
        
        return Tlist

    def gen_Ilist(self):
        Ilist=[] #csr_matrix(eye(self.OrbitalDimension),dtype=complex)]
        
        for d in range(0,self.Dimension):
            D=self.LatticeDimension[d]
            Ilist.append(csr_matrix(eye(D),dtype=complex))
        return Ilist
    
def gen_OnSitePotential(Potential):
    CheckLattice()
    
    """ 
    Generates a Hamiltonian corresponding to on-site potential
    
    The on-site potential A must be an ndarray, where A[nx,ny,..] indicates the on-site potential on site [nx,ny,...]
    """
    
    # Checking for right input"
    
    if not type (Potential)==ndarray or shape(Potential)!=tuple(Lattice.LatticeDimension):
        print(shape(Potential))
        print(Lattice.LatticeDimension)
        raise TypeError("Shape of potential must be ndarray and match physical Hamiltonian")
        
    Vec = ravel(Potential,order="F")
    
    Data=Vec
    Rows=arange(0,Lattice.UnitCells)
    Cols=arange(0,Lattice.UnitCells)
    LatticeMatrix = csr_matrix((Data,(Rows,Cols)),dtype=complex)

    I = csr_matrix(eye(Lattice.OrbitalDimension),dtype=complex)

    return kron(LatticeMatrix,I,format="csr")    
    

def get_density(Vec):
    CheckLattice()
    
    Rho=abs(Vec)**2
    
    Shape=tuple([OrbitalDimension]+list(Lattice.LatticeDimension)+[shape(Rho)[1]])
    RhoMat=reshape(Rho,Shape,order="F")
    
    RhoLattice = sum(RhoMat,axis=0)
    return RhoLattice

    
    
def VecToMatrix(Vec):
    Shape=tuple([OrbitalDimension]+list(Lattice.LatticeDimension))
    
    return reshape(Vec,Shape,order="F")

def MatrixToVec(Mat):
    return ravel(Mat,order="F")    

def get_Identity(format="csr"):
    
   return eye(Hdimension,format=format,dtype=complex)
    
def get_coords(index):
    
    orbital = index%OrbitalDimension
    latticeindex=index//OrbitalDimension
    
    coordlist=[]
    for d in range(0,Dimension):
        coord=(latticeindex//prod(LatticeSize[:d]))%LatticeSize[d]
        coord=coord.astype(int)
        
        coordlist.append(coord)
        
    return [coordlist,orbital]

def CoordinateVec(Dim):
    if Dim + 1 > Dimension:
        raise IndexError("Direction argument must match the lattice dimension")
        
    LatticeIndices = arange(0,Hdimension)//OrbitalDimension
    Vec=(LatticeIndices//prod(LatticeDimension[:Dim]))%LatticeDimension[Dim]
    
    return Vec




