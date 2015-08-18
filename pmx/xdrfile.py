#  -*- mode: python; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*-
# 
#  $Id$
# 
#  Copyright (c) Erik Lindahl, David van der Spoel 2003-2007.
#  Coordinate compression (c) by Frans van Hoesel.
#  Python wrapper (c) by Roland Schulz
# 
#  IN contrast to the rest of Gromacs, XDRFILE is distributed under the
#  BSD license, so you can use it any way you wish, including closed source:
# 
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
# 

from ctypes import *
import os.path

mTrr,mNumPy=1,2
try:
    from numpy import *
    from numpy.ctypeslib import ndpointer
    auto_mode=mNumPy                 
except:
    auto_mode=0


class frame:
    #variables
    #x: rvec*natoms / numpy array if installed
    #box DIM*DIM
    #step 
    #time 
    #prec 
    #lam: lambda

    def __init__(self,n,mode):
        #create vector for x
        if mode&mNumPy:
            self.x=empty((n,3),dtype=float32)
            self.box = empty((3,3),float32)
        else:
            self.x=((c_float*3)*n)() 
            self.box = (c_float*3*3)()


class xdrfile:
    exdrOK, exdrHEADER, exdrSTRING, exdrDOUBLE, exdrINT, exdrFLOAT, exdrUINT, exdr3DX, exdrCLOSE, exdrMAGIC, exdrNOMEM, exdrENDOFFILE, exdrNR = range(13)

    #
    def __init__(self,fn,mode="Auto",ft="Auto"):
        if mode=="NumPy":
          self.mode=mNumPy
          try:
            empty
          except NameError:
              raise IOError("NumPy selected but not correctly installed")
        elif mode=="Std":
          self.mode=0
        elif mode=="Auto":
          self.mode=auto_mode
        else: 
          raise IOError("unsupported mode")
          
        if ft=="Auto":
          ft = os.path.splitext(fn)[1][1:]
          
        if ft=="trr":
            self.mode|=mTrr
        elif ft=="xtc":
            pass
        else:
            raise IOError("Only xtc and trr supported")
        
        #load libxdrfil
        try: 
          p = os.path.join(os.path.dirname(__file__),'_xdrio.so')
          self.xdr=cdll.LoadLibrary(p)
        except:
          raise IOError("_xdrio.so can't be loaded")
          
        #open file
        self.xd = self.xdr.xdrfile_open(fn,"r")
        if not self.xd: raise IOError("Cannot open file: '%s'"%fn)
        
        #read natoms
        natoms=c_int()
        if self.mode&mTrr:
            r=self.xdr.read_trr_natoms(fn,byref(natoms))
        else:
            r=self.xdr.read_xtc_natoms(fn,byref(natoms))
        if r!=self.exdrOK: raise IOError("Error reading: '%s'"%fn)
        self.natoms=natoms.value
        
        #for NumPy define argtypes - ndpointer is not automatically converted to POINTER(c_float)
        #alternative of ctypes.data_as(POINTER(c_float)) requires two version for numpy and c_float array
        if self.mode&mNumPy:
            self.xdr.read_xtc.argtypes=[c_int,c_int,POINTER(c_int),POINTER(c_float),
              ndpointer(ndim=2,dtype=float32),ndpointer(ndim=2,dtype=float32),POINTER(c_float)]
            self.xdr.read_trr.argtypes=[c_int,c_int,POINTER(c_int),POINTER(c_float),POINTER(c_float),
              ndpointer(ndim=2,dtype=float32),ndpointer(ndim=2,dtype=float32),
              POINTER(c_float),POINTER(c_float)]
             
        
    def __iter__(self):
        f = frame(self.natoms,self.mode)
        #temporary c_type variables (frame variables are python type)
        step = c_int()
        time = c_float()
        prec = c_float()
        lam = c_float()
        while 1:
            #read next frame
            if not self.mode&mTrr:
                result = self.xdr.read_xtc(self.xd,self.natoms,byref(step),byref(time),f.box,
                        f.x,byref(prec))
                f.prec=prec.value
            else:
                result = self.xdr.read_trr(self.xd,self.natoms,byref(step),byref(time),byref(lam),
                        f.box,f.x,None,None) #TODO: make v,f possible
                f.lam=lam.value
                
            #check return value
            if result==self.exdrENDOFFILE: break
            if result==self.exdrINT and self.mode&mTrr: 
              break  #TODO: dirty hack. read_trr return exdrINT not exdrENDOFFILE
            if result!=self.exdrOK: raise IOError("Error reading xdr file")
            
            #convert c_type to python 
            f.step=step.value
            f.time=time.value
            yield f
        
