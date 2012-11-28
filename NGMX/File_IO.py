import re, os
from numpy import *
#from numgrompy import xtc_open, xtc_close, xtc_read_headder, xtc_read_coords, xtc_read_coords2, xtc_NT
from numpy.linalg import norm, det, inv
import NGMX_XTC



def _correct_Digits(A, digits):
    """corrects the counting for limited digits in PDB or Gro file:
        given array A
        and the number of digits
        it will look for jumps in numbering of 10**digits and correct after every jump"""
    threshold = 10**digits-1
    A = array(A, dtype='i')
    jumps = nonzero( A[:-1] == A[1:]+threshold )[0] # see positions where count rolls over
    for j in jumps:
        A[j+1:] += 10**digits
    return A

def _LenAng2BOX(Length,Angles):
    la, lb, lc    = Length
    ca, cb, cg    = cos(Angles   *pi/180)
    sg              = sin(Angles[2]*pi/180)
    B =  array([
                [la, lb*cg, lc*cb],
                [0,  lb*sg, lc*(ca-cb*cg)/sg],
                [0,  0,     sqrt(  ( lc )**2 -  ( lc*cb )**2 - ( lc*(ca-cb*cg)/sg )**2  ) ]
                ])
    B[ B < 1e-5] = 0.0 # any component smaller than 10^-6 nm we will consider zero 
    return B

def _BOX2LenAng(BOX):
    Length = sqrt((BOX**2).sum(1))
    B0N = BOX[0]/norm(BOX[0])
    B1N = BOX[1]/norm(BOX[1])
    B2N = BOX[2]/norm(BOX[2])
    Angles = array([
                     arccos(  vdot(B1N, B2N)  ) *180/pi,
                     arccos(  vdot(B0N, B2N)  ) *180/pi,
                     arccos(  vdot(B0N, B1N)  ) *180/pi
                   ])        
    return Length, Angles


class Datafile():
    # protrotype for all files that has all required functions but gives out an error message if accesed
    def __init__(self, filename):
        self.N = self.T = 0 #number of atoms, number of timesteps
        self.filename = filename
        if os.path.isfile(filename):
            self.exsists = True
            try:            
                self._NT_info()
            except:
                raise AssertionError("File corrupted can't read the number of frames or number of atoms")
        else:
            self.exists=False

    def __repr__(self):
        return "<File %s: %d atoms, %d frames>"%(self.filename, self.N, self.T)
        
    def _parse_index(self, index):
        "makes an ndarray index of various kinds of inputs such as strings, tupels, filenames ..."
        from DataObj import Darray
        if type(index) == str:
            if index == 'all':
                index = arange(self.N, dtype=uint64)
            elif os.path.splitext(index)[-1] == 'ndx':
                #parse ndx file
                pass            
            else:
                if len(index.split(',')) != 0:
                    index = index.split(',')
                try:
                    index = array(index, dtype=uint64)
                except:
                    IndexError("String: %s not usable for indexing"%index)                                
        elif type(index) == ndarray or type(index) == Darray:
            if index.dtype is not uint64:
                index = array(index, dtype=uint64)
        
        elif type(index) == list or type(index) == tuple:
            try:
                index = array(index, dtype=uint64)
            except:
                raise IndexError("list or tuple: %s not usable for indexing"%index)
        else:
            raise IndexError("Not supported for indexing")        
        if index.min() < 0 or index.max() > self.N-1:
            raise IndexError("index contains negative numbers or is to large for Data")        
        return index

    def _check_start_stop_skip(self, start, stop, skip):
        try:
            start = int(start);    stop   = int(stop);    skip  = int(skip)
        except:
            raise ValueError("start stop and skip must be integer values: indicating the number of frames in a trajectory (not the time)")
        if abs(start) > self.T:
            raise ValueError("the start value (%d) must be in the range of valid frames (%d)"%(start, self.T))      
        if abs(stop) > self.T:
            raise ValueError("the stop value (%d) must be in the range of valid frames (%d)"%(start, self.T))      
        if skip > self.T:
            raise ValueError("the skip value (%d) must be in the range of valid frames (%d)"%(start, self.T))
        if start  < 0: start  += self.T+1 ## make it possible to index from behind (e.g. stop= -5 is the fifth last frame)
        if stop   < 0: stop   += self.T+1        
        if stop   < start:
            raise ValueError( "the start value (%d) must be larger or equal to the stop value (%d) (Number of Frames = %d)"%(start, stop, self.T) )
        if skip > (stop - start):
            raise ValueError( "the skip frames value %d can't be larger than the difference between start (%s) and stop (%d)"%(skip, start, stop,) )      
        T = int( (stop - start) / skip )      
        return T, start, stop, skip

    def read_structure(self, DataObject=None, index='all'):
        raise IOError("%s: This file does not contain structure information"%self.filename)

    def read_trajectory(self, DataObject=None, index='all', start=1, stop=-1, skip=1):
        raise IOError("%s: This file does not contain trajectory information"%self.filename)

    def iterate_trajectory(self, DataObject=None, index='all', start=1, stop=-1, skip=1):
        raise IOError("%s: This file does not contain trajectory information"%self.filename)

    def write_selection(self, SelObj, start=0, stop=-1, skip=1, Title="numpy Structure", verbose=True):
        raise IOError("%s: writing to this file is not supported yet")

class NDX_File(Datafile):
    def write_selection(self, SelObj, start=0, stop=-1, skip=1, Title="numpy Structure", verbose=True):
        pass

class PDB_File(Datafile):
    #TODO correct unity handling pdb --> Aangstom also for box
    def __init__(self, filename):
        self.reModel = re.compile("^MODEL", re.MULTILINE)
        self.reAtom = re.compile("^(ATOM  |HETATM)(.{5}).(.{4})(.)(.{4})(.)(.{5}).{3}(.{8})(.{8})(.{8})(.{6})(.{6})", re.MULTILINE)
        self.reAcrd = re.compile("^(ATOM  |HETATM).{24}(.{8})(.{8})(.{8})", re.MULTILINE)
        self.reBox  = re.compile( r"""^CRYST1                          #Keyword
                             \s+([\-\.\d]+)  \s+([\-\.\d]+)  \s+([\-\.\d]+) # three box vector length
                             \s+([\-\.\d]+)  \s+([\-\.\d]+)  \s+([\-\.\d]+) # three angles between the vectors
                            """, re.MULTILINE + re.VERBOSE)
        Datafile.__init__(self, filename)
    
    def _NT_info(self):
            s = "\n" + open(self.filename, 'r' ).read()    
            T = len( self.reModel.findall(s) )    # see how many MODEL tags there are
            if T != 0: # if PDB has tags for models                        
                FirstFrame = s.split("\nMODEL")[1]
                N = len( self.reAtom.findall(FirstFrame) )

            else:  #PDB has no Tags for models
                T = 1
                N = len( self.reAtom.findall(s) )
            self.N = N
            self.T = T
            return N,T

    def read_structure(self, DatObj=None, index='all'): #PDB
        from DataObj import Darray
        index = self._parse_index(index)
        s = "\n"+open(self.filename).read()
        T = len(self.reModel.findall(s))
        try:
            if T == 0 or T == 1: #assume only one frame with or without MODEL Tag
                PDB_Data = array( self.reAtom.findall(s) )        
            else:
                PDB_Data = array( self.reAtom.findall(s.split("\nMODEL")[1] ) )
        except:
            raise IOError("Unable to read PDB file maybe file is corrupt..")

        try:
            Length, Angles = array(self.reBox.findall(s)[0], dtype=double).reshape(2,3)
            Box_Data = _LenAng2BOX( Length, Angles )
        except:
            Box_Data = zeros((3,3))
            print("No box info found or corrupt!")
            #raise Warning("No box info found or corrupt!")
            
        ## Data fileds from PDB
        DatObj.BOX       = Darray( Box_Data[newaxis,:,:],                         dtype='f',   dshape = ('T',3,3)   )
        DatObj.RACE      = Darray( PDB_Data[:, 0][index],                         dtype='S6',  dshape = ('N',)      )
        DatObj.ATOMNAME  = Darray( char.strip(PDB_Data[: , 2])[index],            dtype='S4',  dshape = ('N',)      )
        DatObj.ALTLOC    = Darray( PDB_Data[:, 3][index],                         dtype='S1',  dshape = ('N',)      )
        DatObj.RESNAME   = Darray( char.strip(PDB_Data[: , 4])[index],            dtype='S4',  dshape = ('N',)      )        
        DatObj.X         = Darray( PDB_Data[newaxis,:,7:10][:,index],             dtype='f',   dshape = ('T','N',3) )
        DatObj.BFAC      = Darray( PDB_Data[:, 10][index],                        dtype='f',   dshape = ('N',)      )
        DatObj.OCC       = Darray( PDB_Data[:, 11][index],                        dtype='f',   dshape = ('N',)      )
        DatObj.I_ATOM    = Darray( _correct_Digits(PDB_Data[:, 1], 5)[index],     dtype='i',   dshape = ('N',)      )
        DatObj.I_RESIDUE = Darray( _correct_Digits(PDB_Data[:, 6], 4)[index],     dtype='i',   dshape = ('N',)      )
        DatObj.I_CHAIN   = Darray( PDB_Data[:, 5][index],                         dtype='S1',  dshape = ('N',)      )
        DatObj.I_NDX     = Darray( array(index),                                  dtype='i',   dshape = ('N',)      )
        DatObj.I_MODEL   = Darray( ones(DatObj.N),                                dtype='i',   dshape = ('N',)      )


    def read_trajectory(self, DatObj=None, index='all', start=1, stop=-1, skip=1, verbose=True):
        #TODO make the (time) borders right
        T, start, stop, skip = self._check_start_stop_skip( start, stop, skip )        
        if self.N != DatObj.N: raise BufferError("Trajectory you are trying to load does not fit the Structure (Struture %d Atoms; Trajectory %d Atoms)"%(DatObj.nN, self.N)) 
        DatObj.set_T_size(T)
        Frames = ( '\n' + open(self.filename).read() ).split('\nMODEL')[1:]
        if verbose: print("reading file %s"%self.filename)
        for i,f in enumerate(Frames[ start-1 : stop-1 : skip ]):
            if verbose: print"\rreading frame %5d of %5d"%(i+1,T),
            Coord = array(self.reAcrd.findall(f))[:,1:] #get rid of the first "ATOM or HETATOM record"            
            DatObj.X[i] = array( Coord, dtype='f' )
            try:
                Length, Angles = array(self.reBox.findall(f)[0], dtype=double).reshape(2,3)
                DatObj.BOX[i]     = _LenAng2BOX( Length, Angles )
            except:
                DatObj.BOX[i] = zeros((3,3))
                print("No box info found or corrupt!")
                #raise Warning("No box info found or corrupt!")

    def write_selection(self, SelObj, start=0, stop=-1, skip=1, Title="numpy Structure", verbose=True):       
        "writes the current selection as PDB file"        
        if verbose: print("writing selection to %s..."%self.filename)
        pdb_atom_format = "%6s%5d %-4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
        pdb_crst_format = "CRYST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2f\n"
        Title = "TITLE    %s\n"%Title
        self.T = SelObj.nT
        T, start, stop, skip = self._check_start_stop_skip( start, stop, skip )
        
        if  SelObj.nT == 1: #if there is only one time frame
            pdbstr = Title
            pdbstr += "MODEL    1\n"
            if det(SelObj.box) != 0.0:
                Length, Angles =_BOX2LenAng(SelObj.box)
                pdbstr += pdb_crst_format%(Length[0], Length[1], Length[2], Angles[0], Angles[1], Angles[2])
            zipper = zip(SelObj.race, SelObj.i_atom, SelObj.atomname, SelObj.altloc, SelObj.resname, SelObj.i_chain,
                        SelObj.i_residue, SelObj.x[:,0], SelObj.x[:,1], SelObj.x[:,2], SelObj.bfac, SelObj.occ )
            for race, iatom, name, altloc, resname, ichain, iresidue, x, y, z, bfac, occ in zipper:
                pdbstr += pdb_atom_format%(race, iatom%10**5, name, altloc, resname, ichain, iresidue%10**4, x, y, z, bfac, occ )
            pdbstr += "ENDMDL\n"
                
        else: # if more than one time frames
            # get the string ready that will be the same for every timeframe (ATOM, resname number..) only coordinates will change
            pdb_atom_format = "%6s%5d %-4s%1s%3s %1s%4d    %s%s%s%6.2f%6.2f\n" #change the %8.3f to %s to insert by % later 
            zipper = zip(SelObj.race, SelObj.i_atom, SelObj.atomname, SelObj.altloc, SelObj.resname, SelObj.i_chain,
                        SelObj.i_residue, SelObj.bfac, SelObj.occ )
            frame_const = ''
            for race, iatom, name, altloc, resname, ichain, iresidue, bfac, occ in zipper:
                frame_const += pdb_atom_format%(race, iatom%10**5, name, altloc, resname, ichain, iresidue%10**4, "%8.3f", "%8.3f", "%8.3f", bfac, occ )
            
            
            pdbstr=''
            # loop over every frame            
            for ii, i in enumerate(range(SelObj.nT)[start:stop:skip]):                
                try: # try to get the timeinfo into the title
                    #Tit_time = Tit_time[:-1]+"\tt=   %8.3f\n"%SelObj.time[i]
                    Tit_time = Title[:-1] + "\tt=   %8.3f\n"%SelObj.time[i]
                except:
                    Tit_time = Title  
                if verbose: print "\r...writing time frame %i of %i"%(ii+1, T),
                pdbstr += Tit_time
                pdbstr += "MODEL    %5i\n"%i
                # write box info
                if det(SelObj.box[i]) != 0.0:
                    Length, Angles =_BOX2LenAng(SelObj.box[i])
                    pdbstr += pdb_crst_format%(Length[0], Length[1], Length[2], Angles[0], Angles[1], Angles[2])

                #write coordiantes
                pdbstr += frame_const%tuple( SelObj.x[i].flatten() )
                #write endmodel to close frame
                pdbstr += "ENDMDL\n"            
        open(self.filename, 'w').write(pdbstr)
        if verbose: print("...done")

    def write_selectionII(self, SelObj, start=0, stop=-1, skip=1, Title="numpy Structure", verbose=True):       
        "writes the current selection as PDB file easier but slower... still really slow maybe port to C..."        
        if verbose: print("writing selection to %s..."%self.filename)
        pdb_atom_format = "%6s%5d %-4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
        pdb_crst_format = "CRYST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2f\n"
        Title = "TITLE    %s\n"%Title
        self.T = SelObj.nT
        T, start, stop, skip = self._check_start_stop_skip( start, stop, skip )
        print "start %i, stop %i, skip %i"%(start,stop, skip)
        fh = open(self.filename, 'w')
        pdbstr = Title
        SelObj = SelObj.T[start:stop:skip]
        for i,frame in enumerate(SelObj.T):
            if verbose: print "\rframe %i of %i"%(i,SelObj.nT),
            pdbstr += "MODEL    %5i\n"%i ## Model
            if det(frame.box) != 0.0:    ## Cryist
                Length, Angles =_BOX2LenAng(SelObj.box[i])
                pdbstr += pdb_crst_format%(Length[0], Length[1], Length[2], Angles[0], Angles[1], Angles[2])
            
            ### atom records    
            pdbdata = zip(frame.race, frame.i_atom, frame.atomname, frame.altloc, frame.resname, frame.i_chain,
                        frame.i_residue, frame.x[:,0], frame.x[:,1], frame.x[:,2], frame.bfac, frame.occ )
            for race, iatom, name, altloc, resname, ichain, iresidue, x, y, z, bfac, occ in pdbdata:
                pdbstr += pdb_atom_format%(race, iatom%10**5, name, altloc, resname, ichain, iresidue%10**4, x, y, z, bfac, occ )
            pdbstr += "ENDMDL\n"
            fh.write(pdbstr)
            pdbstr = ''             
        fh.close()
        if verbose: print("...done")
        
        
    def append_selection(self, SelObj, start=0, stop=-1, skip=1, new_Model=False, verbose=True):
        """ append an object to an existing pdb File"""
        
        if verbose: print("appending selection to %s..."%self.filename)
        pdb_atom_format = "%6s%5d %-4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
        T, start, stop, skip = self._check_start_stop_skip( start, stop, skip )        
        Fstr = ('\n' + open(self.filemane,'r').read()).split("\nENDMDL")
                
        if  SelObj.T == 1: #if there is only one time frame           
            
            
            zipper = zip(SelObj.race, SelObj.i_atom, SelObj.atomname, SelObj.altloc, SelObj.resname, SelObj.i_chain,
                        SelObj.i_residue, SelObj.x[:,0], SelObj.x[:,1], SelObj.x[:,2], SelObj.bfac, SelObj.occ )
            for race, iatom, name, altloc, resname, ichain, iresidue, x, y, z, bfac, occ in zipper:
                pdbstr += pdb_atom_format%(race, iatom%10**5, name, altloc, resname, ichain, iresidue%10**4, x, y, z, bfac, occ )
            pdbstr += "ENDMDL\n"
                
        else: # if more than one time frames
            # get the string ready that will be the same for every timeframe (ATOM, resname number..) only coordinates will change
            pdb_atom_format = "%6s%5d %-4s%1s%3s %1s%4d    %s%s%s%6.2f%6.2f\n" #change the %8.3f to %s to insert by % later 
            zipper = zip(SelObj.race, SelObj.i_atom, SelObj.atomname, SelObj.altloc, SelObj.resname, SelObj.i_chain,
                        SelObj.i_residue, SelObj.bfac, SelObj.occ )
            frame_const = ''
            for race, iatom, name, altloc, resname, ichain, iresidue, bfac, occ in zipper:
                frame_const += pdb_atom_format%(race, iatom%10**5, name, altloc, resname, ichain, iresidue%10**4, "%8.3f", "%8.3f", "%8.3f", bfac, occ )
            # loop over every frame
            for i in range(SelObj.T)[start:stop:skip]:
                if verbose: print( "\r...writing time frame %i of %i"%(i+1, T) )
                pdbstr += "MODEL    %5i\n"%i
                # write box info
                if det(SelObj.box[i]) != 0.0:
                    Length, Angles =_BOX2LenAng(SelObj.box[i])
                    pdbstr += pdb_crst_format%(Length[0], Length[1], Length[2], Angles[0], Angles[1], Angles[2])
                #write time info
                try:
                    pdbstr += "t=   %8.3f\n"%SelObj.time[i]
                except:
                    pass # if there is no time info ... don't put it in the pdb file
                #write coordiantes
                pdbstr += frame_const%tuple( SelObj.x[:,i].flatten() )
                #write endmodel to close frame
                pdbstr += "ENDMDL\n"            
        open(self.filename, 'w').write(pdbstr)
        if verbose: print("...done")

class GRO_File(Datafile):
    def __init__(self, filename):
        self.reBox_triclinic = re.compile("^\s+([\.-0-9]+)\s+([\.-0-9]+)\s+([\.-0-9]+)\s*\n", re.MULTILINE)
        self.reBox_non_tric  = re.compile("^\s+([\.-0-9]+)\s+([\.-0-9]+)\s+([\.-0-9]+)\s+([\.-0-9]+)\s+([\.-0-9]+)\s+([\.-0-9]+)\s+([\.-0-9]+)\s+([\.-0-9]+)\s+([\.-0-9]+)\s*\n", re.MULTILINE)
        self.reAtom = re.compile("^(ATOM  |HETATM)(.{5}).(.{4})(.)(.{4})(.)(.{5}).{3}(.{8})(.{8})(.{8})(.{6})(.{6})", re.MULTILINE)
        self.reAcrd = re.compile("^(ATOM  |HETATM).{5}..{4}..{4}..{5}.{3}(.{8})(.{8})(.{8}).{6}.{6}", re.MULTILINE)
        Datafile.__init__(self, filename)
    def _NT_info(self):
        pass

class XTC_File(Datafile):
    #TODO if file does not exist --> errormessage insted of segfalut!
    def _NT_info(self):
        self.N, self.T = NGMX_XTC.xdrfile_NT(self.filename) #C function in NGMX_XTC->xdrfile_NT
        return self.N, self.T    

    def read_trajectory(self, DatObj=None, index='all', start=0, stop=-1, skip=1, verbose=True, unity='A'):
        "reads a xtc and writes the coordinates, box info, time and step to the appropriate arrays"
        #make sure that T and N are set right otherwise we could run out of the memory space        
        T, start, stop, skip = self._check_start_stop_skip( start, stop, skip )
        index = self._parse_index(index)
            
        #for long files it may be good to see that something is happening...
        if verbose:    verbose = 1
        else:        verbose = 0
        
        #numbers have to be multiplied anyway so it is quicker to directly do that in the C code
        if      unity=='A': scaling  = 10.0
        elif    unity=='nm': scaling =  1.0
        elif    (type(unity) == int) or (type(unity)) == float: scaling = float(scaling)
        else:   raise TypeError("Only supports 'A' (Aangstroms) and 'nm' (nano meters) or a number as the scaling factor with respect to nm.")
        
        if not os.path.isfile(self.filename):
            raise AssertionError("File %s does not exist!"%self.filename)
            
        DatObj.set_T_size(T)
        NGMX_XTC.xdrfile_read(self.filename, DatObj.X, DatObj.BOX, DatObj.TIME, DatObj.STEP, index, start, stop, skip, scaling, verbose)            
        DatObj.unity = 'nm'
            
    def write_selection(self, DatObj=None, start=0, stop=-1, skip=1, verbose=True, unity='A'):
        #make sure that T and N are set right otherwise we could run out of the memory space
        self.N = DatObj.nN; self.T = DatObj.nT
        T, start, stop, skip = self._check_start_stop_skip( start, stop, skip )
        D    = DatObj.T[start:stop:skip] #make a dataobjet, that only carries the data we want to write out
        self.T = T; self.N = D.nN
        #T, start, stop, skip = self._check_start_stop_skip( start, stop, skip )        
        #for long files it may be good to see that something is happening...
        if verbose:    verbose = 1
        else:        verbose = 0        
        
        if      unity=='A': scaling  = 10.0
        elif    unity=='nm': scaling =  1.0
        elif    (type(unity) == int) or (type(unity)) == float: scaling = float(scaling)
        else:    raise TypeError("Only supports 'A' (Aangstroms) and 'nm' (nano meters) or a number as the scaling factor with respect to nm.")             
        NGMX_XTC.xdrfile_write(self.filename, D.x, D.box, D.time, D.step, scaling, verbose)
        


KnownFiles = { '.pdb' : PDB_File,
               '.xtc' : XTC_File
             }
