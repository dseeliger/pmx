from numpy import *
import os, re
import library as lib
from DataObj import Darray



## global functions
def unique_retord(A):
    "unique but sorted!"
    try:
        Ad, Ai = unique(A, return_index=True)
    except: 
        #for some wierd reason numpy unique just stopped to work with the pointer so we use a has if this fails here
        Ad, Ai = unique( array([ hash(a) for a in A ]), return_index=True)    
    Ai.sort() #retain order!        
    return A[ Ai ]


class N_Slicer():
    """ The N_slicer is a class used to slice the Data along all Axes that sclae with the number of Atoms in the system.
Here we can address arbitrary positions in the DataObject by giving an array that would index """
    l  = property(lambda self: len(self), doc="How many along this axis")
    def __init__(self, SelObj, crt):        
        self.SelObj = SelObj
        self.scale = 'N'  
        if type(crt) is ndarray or type(crt) is Darray:
            self.crt = crt
            #This was implemented to have them in order, but maybe we don't want that all the time...
            #self.crt = unique(crt) # let's see how this affects performance....
            #self.crt.sort()
             
        elif type(crt) is slice:
            self.crt = arange(self.SelObj.D.N)[crt] #convert slices to an indexing array
            
        elif type(crt) is int:
            self.crt = array([crt])
            
        elif crt == 'all':        #convinience if all Atoms should be selected 
            self.crt =  arange(self.SelObj.D.N)
        else:
            raise TypeError("Indexing array for N is not valid!")
 
    def __len__(self):
        return len(self.crt)
    
    def __getitem__(self, slc):
        newN = self.crt[slc]
        newN = array([newN]).flatten()  # to cope with integers if single element is selected      
        return SelObj( self.SelObj.D, newN, self.SelObj.T.crt )

    def __iter__(self):
        SO = SelObj(self.SelObj.D, 0, self.SelObj.T.crt)              
        for i in self.crt: #for every time frame number
            SO.N.crt = i
            yield SO        
    
    def is_empty(self):
        return len(self) == 0

class T_Slicer():
    "Here we can slice the T dimension with a sclice. if an indexing array is given we will try to convert it to a slice"
    l  = property(lambda self: len(self) , doc="How many along this axis")    
    def __init__(self, SelObj, crt):
        self.SelObj = SelObj
        self.scale = 'T'        
        if type(crt) is slice: # this is how it should be a slice for T
            start, stop, step = crt.indices(self.SelObj.D.T)
            if crt.start is None: start = None              #if it was None keep it that way!
            elif crt.start < 0: start -=  self.SelObj.D.T   #if it was negative keep it that way!            
            if crt.stop  is None: stop  =  None
            elif crt.stop  < 0: stop  -=  self.SelObj.D.T
            self.crt = slice(start, stop, step)                        
        elif type(crt) is int:
            self.crt = crt
        elif crt == 'all':          # convinience if all time frames should be selected 
            self.crt =  slice(None)
        elif type(crt) is ndarray: # if someone tried to use an indexing array we will try to recover it... 
            if len(crt) == 1:       # if it is only one number e.g. just one time frame we will just take that
                self.crt = crt[0]
            else:                   # we will try to make a valid slice out of it...
                self.crt = self.ndxarray_2_sclice(crt)
                
    def __len__(self):
        if  type(self.crt) is int:
            return 1
        else:
            return self.slice_len(self.crt ) 
 
    def slice_len(self, slc):
        start, stop, step = slc.indices(self.SelObj.D.T)
        return int(ceil(  (stop-start) / float(abs(step))  ))
            
    def ndxarray_2_sclice(self, ndxarray):
        "tries to convert a sequence of index into a slice object and gives out an error if this fails"
        ndxarray = unique(ndxarray); ndxarray.sort()
        start = ndxarray[0]; stop = ndxarray[-1]+1; step = ndxarray[1] - ndxarray[0]
        A = arange(start, stop, step) == ndxarray   ## if this test holds the array can be represented as a slice
        if type(A) is ndarray: A = A.all()     #if it was two arrays of same length we will get a boolean vector
        if A : 
            return slice(start, stop, step)
        else:
            raise TypeError("Indexing array for T can not be converted into valid slice(start, stop, step)! You can only index T in a regular way! ")
 
    def __getitem__(self, slc):
        # if indexed by a number
        if type(slc) is int:
            start, stop, step = self.crt.indices(self.SelObj.D.T)
            i = start + slc*step
            return SelObj( self.SelObj.D, self.SelObj.N.crt, slice(i,i+1,None) )
        
        elif type(slc) is slice :
            ### combine the two slicing operators into one
            Start1, Stop1, Step1 = self.crt.indices(self.SelObj.D.T)
            Start2, Stop2, Step2 = slc.indices(len(self))
            
            # get the new step
            newStep = Step1 * Step2
             
            # get the new start               
            if Start1 is None and Start2 is None: newStart = None
            else:
                if Start1 is None:                    newStart = Start2
                elif Start2 is None:                  newStart = Start1                
                else:                                 newStart = Start1 + (Start2 * Step1)
                if self.crt.start < 0:                newStart - self.SelObj.D.T
                
            #get the new end
            if Stop1 is None and Start2 is None: newStart = None
            else:
                if Stop1 is None:                     newStop = Stop2
                elif Start2 is None:                  newStop = Stop1                
                else:                                 newStop = Start1 + (Stop2 * Step1)
                if self.crt.stop < 0:                 newStop - self.SelObj.D.T
                
            newSlc = slice(newStart, newStop, newStep)
            return SelObj( self.SelObj.D, self.SelObj.N.crt, newSlc )
        
        else:
            raise IndexError("Can only slice the Time dimension ('[start:stop:skip]') or give the number of a specific time frame within selection")

    def __iter__(self):
        SO = SelObj(self.SelObj.D, self.SelObj.N.crt, 0)              
        for i in arange( self.SelObj.D.T )[self.crt]: #for every time frame number
            SO.T.crt = slice(i,i+1)
            yield SO
    
    def is_empty(self):
        return len(self) == 0


"""ideas: random_residue(), move_to(), cartesian_box(), closest_to(xyz)"""

class SelObj(object):
    nN      = property(lambda self: self.N.l, doc="How many atoms in selection")
    nT      = property(lambda self: self.T.l, doc="How many time frames in selection")    
    xNaxis  = property(lambda self: 0 if self.nT == 1 else 1) #property to work with the squeezed self.x array (returns the axis of N ---> mean(self.xNaxis) mean along the N axis)
    xCaxis  = property(lambda self: 1 if self.nT == 1 else 2) #property to work with the squeezed self.x array (returns the axis of N ---> mean(self.xNaxis) mean along the N axis)

    
    
    def __init__(self, D, Ncrt, Tcrt=slice(None,None,None), verbose=True):
        from DataObj import Darray, Oarray
        self.D = D                                  # object where all Data is stored
        self.N = N_Slicer(self, Ncrt)            # selection of atoms in D
        self.T = T_Slicer(self, Tcrt)            # selection of time frames in D
        if self.N.is_empty() or self.T.is_empty():
            raise MemoryError("Error, will not create empty Selection")
        
        for attr, [obj, dim, init, dat_type]  in self.D.attributes.items():
            if obj is Darray:
                setattr(self.__class__, attr.lower(),       property( **self._Darray_ppt(attr) ))
                setattr(self.__class__, attr.lower()+"_uf", property( **self._Darray_ppt_unflat(attr) ))
            elif obj is Oarray:
                setattr(self.__class__, attr.lower(), property( **self._Oarray_ppt(attr) ) )            

    def __slice(self, dim):
        "converts a tuple with dimensions like ('N','T',3) into a tuple to slice the self.D.attr Array; it implicitly flattens the array"            
        rdim=[]        
        for d in dim:
            if d == 'N':
                if self.nN == 1: #if we have only one Atom flatten the data so that assignments will fit (no hidden dimensions)
                    rdim.append(self.N.crt[0])
                else:
                    rdim.append( self.N.crt )
            elif d == 'T':
                if self.nT == 1: #if we have only one time frame select that time frame with an integer so that assignments will fit (no hidden dimensions) --> implicid squeeze
                    if self.T.crt.start is None:
                        rdim.append( 0 )
                    else:
                        rdim.append( self.T.crt.start )
                else:                    
                    rdim.append( self.T.crt )
            elif type(d) is int: #if the dimension is an integer (like 3 for three dimensions) select the whole range with a "none" slice
                if d == 1:
                    rdim.append( 0 )
                else:
                    rdim.append( slice(None) )
            else:
                raise TypeError(  "Error in attribute: dimensionality may only contain 'N','T' or a fixed number but is %s"%(str(dim) )  )
        return tuple(rdim)

    def __slice_unflat(self, dim):
        "converts a tuple with dimensions like ('N','T',3) into a tuple to slice the self.D.attr Array; it will slice an unflattened array"
        rdim=[]       
        for d in dim:
            if d == 'N':
                rdim.append( self.N.crt )
            elif d == 'T':
                rdim.append( self.T.crt )
            elif type(d) is int: #if the dimension is an integer (like 3 for three dimensions) select the whole range with a "none" slice
                rdim.append( slice(None) )
            else:
                raise TypeError(  "Error in attribute: dimensionality may only contain 'N','T' or a fixed number but is %s"%(str(dim) )  )
        return tuple(rdim)

    def _onN(self):
        """returns a slice object, that can be used on arrays to fill up the axis, so that the array will align with N axis"""
        if self.nT != 1:
            return (newaxis, slice(None)) # [:]
        else:
            return slice(None) # [newaxis,:]
    onN     = property( _onN )

    def _inN(self):
        """returns a slice object, that will insert an N dimension to subsitute the N dimension of the coordinates array if it exists"""
        if self.nT != 1:
            return (slice(None),)


    def _axN(self):
        if not self.nT == 1:
            return 1 
        else:
            return 0
    axN     = property( _axN )

    def _onC(self):
        """returns a slice object, that can be used on arrays to fill up the T or N axis if they exsits"""
        if self.nT != 1 and self.nN != 1:           #both are not one
            return (newaxis, newaxis, slice(None))
        elif (self.nT != 1) != (self.nN != 1):      # one of them is one
            return (newaxis, slice(None))
        else:                                       #both are one
            return slice(None)
    onC     = property( _onC )        

    def _axC(self):
        if self.nT != 1 and self.nN != 1:           #both are not one
            return 2
        elif (self.nT != 1) != (self.nN != 1):      # one of them is one
            return 1
        else:                                       #both are one
            return 0
    axC     = property( _axC )


                
    def _Darray_ppt(self, atr_name):   ### Contains Data like   
        doc = "pointer to all the array lines for %s for this selection"%atr_name        
        def fget(self):
            DA = getattr( self.D, atr_name.upper() )
            return DA[ self.__slice(DA.dshape) ]# .squeeze()  # squeeze is implicid in __slice()
        def fset(self, value):
            # TODO -> for T=1 try unflatten the value to fit Darray!
            DA = getattr( self.D,atr_name.upper() )            
            DA[ self.__slice(DA.dshape) ] = value
        return {'fget':fget, 'fset':fset, 'doc':doc}

    def _Darray_ppt_unflat(self, atr_name):
        doc = "pointer to all the array lines for %s for this selection, keeping all dimensions even empty ones"%atr_name
        def fget(self):
            DA = getattr( self.D, atr_name.upper() )
            return DA[ self.__slice_unflat(DA.dshape) ]# .squeeze()  # squeeze is implicid in __slice()
        def fset(self, value):
            # TODO -> for T=1 try unflatten the value to fit Darray!
            DA = getattr( self.D,atr_name.upper() )            
            DA[ self.__slice_unflat(DA.dshape) ] = value
        return {'fget':fget, 'fset':fset, 'doc':doc}


    def _Oarray_ppt(self, atr_name):   ### Object properties, that help navigate through the Model
        doc = "pointer to all the array lines for %s for this selection"%atr_name        
        def fget(self):
            OA = getattr( self.D, atr_name.upper() )
            R = unique_retord( OA[ self.__slice(OA.dshape) ] ) #only return one pointer to every object
            if len(R) == 1: R = R[0] #if there is only one object just return that
            return R
        def fset(self, value):
            # maybe check if value is a selection object ???
            DA = getattr( self.D,atr_name.upper() )
            DA[ self.__slice(DA.dshape) ] = value
        return {'fget':fget, 'fset':fset, 'doc':doc}
    
    def _D_test(self, selection):
        if self.D is not selection.D:
            raise MemoryError("Operation can't be done: the two objects don't share the dame data object")

    def __repr__(self, O_type="Selection"):
        #c = r = 0; a = 1
        ret_str = "<%s:"%O_type    
        try:
            if self.chain.__class__.__name__ == 'Chain': # can't check for type --> circular imports ...
                ret_str += " of Chain '%s';"% self.chain.i_chain[0]
            else:
                ret_str += " %d Chains;"% len(self.chain)
        except:
            pass
        try:
            if self.residue.__class__.__name__ == 'Residue':
                ret_str += " of Residue %s %s;"% (self.residue.resname[0], self.residue.i_residue[0])
            else:
                ret_str += " %d Residues;"% len(self.residue)
        except:
            pass
        try:
            if self.i_atom.shape == ():
                ret_str += " single Atom '%s';"%self.atomname
            else:
                ret_str += " %d Atoms;"% len(self.i_atom)
        except:
            pass
        if self.nT != 1:            
            ret_str = ret_str[:-1] + ' | %d frames '%self.nT        
        return ret_str[:-1] +'>'
          
    def __len__(self):
        return (self.N.l, self.T.l)     #return a tuple of (Number of Atoms, Number of time frames)

    def __and__(self, B):
        self._D_test(B)
        newN_crt = self.N.crt[ in1d(self.N.crt, B.N.crt) ]; newN_crt.sort()        
        return SelObj(  self.D, newN_crt, self.T.crt  )

    def __mult__(self,B):
        self._D_test(B)
        newN_crt = self.N.crt[ in1d(self.N.crt, B.N.crt) ]; newN_crt.sort()        
        return SelObj(  self.D, newN_crt, self.T.crt  )        
    
    def __or__(self, B):
        self._D_test(B)
        newN_crt = unique( r_[self.N.crt, B.N.crt]); newN_crt.sort()
        return SelObj(  self.D, newN_crt, self.T.crt  )

    def __add__(self, B):
        self._D_test(B)
        newN_crt = unique( r_[self.N.crt, B.N.crt]); newN_crt.sort()
        return SelObj(  self.D, newN_crt, self.T.crt  )

    def __sub__(self, B):
        self._D_test(B)
        newN_crt = self.N.crt[ -in1d(self.N.crt, B.N.crt) ]; newN_crt.sort()
        return SelObj(  self.D, newN_crt, self.T.crt  )

    def __contains__(self, B):
        self._D_test(B)
        return in1d(B.N.crt, self.N.crt).all()
    
    def __neg__(self):
        AA = arange(self.D.N) # All atoms
        return SelObj(  self.D, setdiff1d(AA, self.N.crt), self.T.crt  )

    def _make_special_selections(self, verbose=False):
        "Tries to make sense of the objects by pre-assigning selections like Protein, virtual sites etc."
        # Protein
        try:
            prot_crt = nonzero( in1d(self.resname, lib._protein_residues) )[0]
            self.protein = SelObj( self.D, prot_crt )
        except MemoryError:
            if verbose: print "No protein found!"
        # Hydrogens
        try:
            self.hydrogens =self.fetch("H", 'by_symbol' )
        except:
            if verbose: print "No hydrogens found!"
        # Virtual sites
        try:
            if hasattr(self, "symbol"):
                vsite_crt = nonzero(self.symbol == 'DV')[0]
            else:
                vsite_crt = nonzero( in1d(self.atomname, lib._virtual_sites) )[0]
            self.virtual_sites = SelObj( self.D, vsite_crt, self.T.crt )
        except MemoryError:
            if verbose: print "No virtual sites found!"        
        #Lipids
        try:
            lipid_crt = nonzero( in1d(self.resname, lib._lipid_residues) )[0]
            self.lipids = SelObj( self.D, lipid_crt, self.T.crt )
        except MemoryError:
            if verbose: print "No lipids found!"            
        #Water
        try:
            water_crt = nonzero( in1d(self.resname, lib._water) )[0]
            self.water = SelObj( self.D, water_crt, self.T.crt  )
        except MemoryError:
            if verbose: print "No water found!"            
        #ions
        try:
            ion_crt = nonzero( in1d(self.resname, lib._ions) )[0]
            self.ions = SelObj( self.D, ion_crt, self.T.crt  )
        except MemoryError:
            if verbose: print "No ions found!"            
        #Nuclic acids
        try:
            nucacid_crt = nonzero( in1d(self.resname, lib._nucleic_acids) )[0]
            self.nucleic_acids = SelObj( self.D, nucacid_crt, self.T.crt  ) 
        except MemoryError:
            if verbose: print "No nuclic acids found!"

#### comparative boolean operations
    def has_same_NAtoms(self, ref):
        return self.nN == ref.nN
        
    def has_same_or_1_Atoms(self, ref):
        return ref.nN in [self.nN, 1]
        
    def has_same_NFrames(self, ref):
        return self.nT == ref.nT

    def has_same_or_1_Frames(self, ref):
        return ref.nT in [self.nT, 1]


#### Search and database functions

    def fetch(self, key=None, how='by_atomname', inv=False, wildcard=False):
        """fetch atoms by selection,
        key: give a list of atoms to find e.g. ['CA','CB',N] ; 'CA'
        how: tell weather you want to search in e.g. 'by_resname'; 'by_atomname'; 'by_symbol'; 'by_race' ...
        inv: invert selection"""
        from DataObj import Darray        
    
        def search_seq(pattern):
            from DataObj import Oarray
            seq = self.sequence()
            seq = re.sub('<.*?>','0',seq) # get rid of all unknown <xxx> residures and replace them with X so that string position matches up with position in sequence 
            crt = set()
            if wildcard: re.sub('X','.',pattern)            
            for Snip in re.finditer(pattern.upper(), seq): #use pattern as regular expression to find all ocuraces of that snipplet in the sequence
                b_slice= slice(None) # to only get the slice in the brackets
                try: # see if there are any brackets in the selection
                    b_op = pattern.index('(')
                    b_cl = pattern.index(')')
                    if b_op > b_cl: raise SyntaxError("If you provide brackets in the search pattern you can only provide one pair and do it right!") #see if they are set right
                    b_slice = slice(b_op,b_cl-1)
                except ValueError:
                    pass # no brackets are geien or in a wrong way ...
                [ crt.update( r.N.crt ) for r in  self.residue[ Snip.start(): Snip.end() ][b_slice] ]
            crt = array( list(crt) ); crt.sort()
            return SelObj( self.D, array(crt), self.T.crt )
        
        def search_in(key, in_array):
            if not (type(key) is Darray or type(key) is list): key = [key] 
            if not wildcard:
                sele = nonzero( in1d( in_array , key) )[0]
            else:
                sele = zeros(  len( self.N )  ) > 0 # vector with all False of len N
                for k in key:
                    sele +=  char.find( in_array , k ) # vector with Ture where key matched        
            return SelObj(self.D, self.N.crt[sele], self.T.crt)        
        #error handling
        try:
            if how[:3] != "by_": raise Exception()
            field = how[3:]
            if not hasattr(self, field): raise LookupError("selection does not have a property %s"%field)                    
        except:
            raise LookupError("the search method 'how' must start with 'by_' followed by the field you want to search in. %s is an invalid selection"%how)
        # see if we search by field or by sequence
        if field == "sequence":
            sel = search_seq(key)
        else:
            try:
                in_array =  getattr(self, field)
                if type(in_array) is not Darray: raise LookupError("you can't search in %s, it is not a valid array"%field)
                sel = search_in(key, in_array)
            except:
                raise LookupError("Unable to find any %s in %s"%(key, field))        
        if inv: 
            return (self - sel) # return the objects not within the selection
        else:
            return sel

    def sequence(self):
        "generate the one letter sequence"
        from DataObj import Oarray
        seq_str = ['']
        if type(self.residue) is not Oarray:        #if there is only one residue
            AA = Oarray( (1,) )
            AA[:] = self.residue
        else:                                       #if there are multiple residues (what should be the case...)
            AA = self.residue       
        for aa in AA:
            try:
                aa_name = aa.resname[0]
                seq_str += lib._one_letter[aa_name]
            except:
                seq_str += '<%s>'%aa_name        
        return ''.join(seq_str)

               
    def dist(self, reference, squared=False, unflat=False, pbc=False):        
        """returns the distance matrix of every atom to every other atom of the reference"""
        if not self.has_same_or_1_Frames(reference): raise AssertionError("reference must have the same number of timeframes or only one timeframe")                        
        D = self.x_uf[:,newaxis] - reference.x_uf[:,:,newaxis]
        
        if pbc: # takes care of distances that are closer over the periodic image
            halfbox = self.box.diagonal(axis1=1, axis2=2)/2 #So far only rectangular boxes...
            D       = abs(D) % halfbox[:,newaxis,newaxis]        
        
        D = (D**2).sum(3) # x**2 + y**2 + z**2
        
        if not squared: # does the expensive sqrt operation if not unselected by the user
            D = sqrt(D)
        
        if not unflat: # gets rid of dimensions that are one...
            D = D.squeeze()
                   
        return D # return the big matrix!        


    def dist_residue(self, reference, squared=False, unflat=False, pbc=False):        
        """returns the distance matrix of every atom to every other atom of the reference"""
        if not self.has_same_or_1_Frames(reference): raise AssertionError("reference must have the same number of timeframes or only one timeframe")                               
        D = self.x_uf[:,newaxis] - reference.x_uf[:,:,newaxis]
        
        if pbc: # takes care of distances that are closer over the periodic image
            halfbox = self.box.diagonal(axis1=1, axis2=2)/2 #So far only rectangular boxes...
            D       = abs(D) % halfbox[:,newaxis,newaxis]        
        
        D = (D**2).sum(3) # x**2 + y**2 + z**2
        
        if not squared: # does the expensive sqrt operation if not unselected by the user
            D = sqrt(D)
        
        if not unflat: # gets rid of dimensions that are one...
            D = D.squeeze()
                   
        return D # return the big matrix!





#    def closest(self, center=None, N_closest=10, farthest=False, ret_distance=False):
#        """ Returns the N closest Atoms to the given center as a selection
#            --> farthest=True selects the farthest residues instead
#            --> ret_distance retruns the distance in an array
#        """            
#        if center.shape not in [ (self.nT,3), (3,) ]:
#            raise AssertionError("there must be one coordinate or one for every time frame!")
#        
#        if self.nT == len(center):
#            D = ((self.x - center[:,newaxis])**2).sum(self.xCaxis)
#        else
#            D = ((self.x - center)**2).sum(self.xCaxis)
#        crt = self.N.crt[argsort(D)[] ]
        
    def by_residue(self):
        "Return a selection of complete residues (if they are in the Model!)"        
        try:
            return self.residue.asSelection().T[ self.T.crt ]
        except AttributeError:
            return self.residue.T[ self.T.crt ]

    def by_chain(self):
        "Return a selection of complete chains"
        try:
            return self.chain.asSelection().T[ self.T.crt ]
        except AttributeError:
            return self.chain.T[ self.T.crt ]

    def newModel(self, verbose=True):
        """Retruns a Model of the same Selection but with a new DataObj ('.D') behind it containing only the atoms in the selection.
            Use this to reduce the number of atoms prior to reading in a trajectory!"""
        from DataObj import DataObj, Oarray, Darray        
        kwarg = {} #get all attributes of the old D object
        Oarray_fields = {}
        for Attr, [obj, dim, init, dat_type] in self.D.attributes.items():
            Attr = Attr.upper()
            if obj is Darray:
                Dat = getattr( self.D, Attr )
                Dat = Dat[ self.__slice_unflat(Dat.dshape) ]
                #if type(Dat) is not Darray: #if it only has one element the result will be that element, so we have to transfer it back to an Darray
                #    Dat = Darray( array([Dat]), dshape=getattr( self.D, Attr ).dshape )                    
                kwarg.update({ Attr: Dat })
            elif obj is Oarray:
                Oarray_fields.update({Attr: [obj, dim, init, dat_type]})
        newD = DataObj(self.nN, self.nT, [Oarray_fields], **kwarg)
        newD.parse_data( verbose )
        M =  Model( DataObj = newD ) 
        M._make_special_selections()  
        newD.MODEL[:] = M
        return M

    def newModel_Slice(self, index, T_slice=None, N_slice=None, verbose=True):
        """Slices the Data in N and T direction according to supplied Slicing indexes
        those can be:
        1)lists/arrays with the indexes that we want from this direction [1,5,2,27...]
        2)an array  of the length of the dimension with True/False
        3)a slicing object slice(10,100,2) to go from the 10th to the 100th in steps of two
        
        ... and will return an newObject (not a selection) with this data"""
        
        from DataObj import DataObj, Oarray, Darray        
        kwarg = {} #get all attributes of the old D object
        Oarray_fields = {}
        for Attr, [obj, dim, init, dat_type] in self.D.attributes.items():
            Attr = Attr.upper()
            if obj is Darray:
                Dat = getattr( self.D, Attr )
                
                # start slicing with the larger dimension
                if self.nN > self.nT: slcOrder = [['N',N_slice],['T',T_slice]]
                else: slcOrder = [['T',T_slice],['N',N_slice]]                               
                # do the slicing on the data
                for d, d_slice in slcOrder:               
                    try:
                        iN = Dat.dshape.index(d)        # find where the datashape is N/T
                        iN = [slice(None)]*iN           # make as many empty slices till we get to the dimension we really want to slice
                        Dat = Dat[ tuple(iN + [d_slice])   ] # slice the array
                    except ValueError:
                        pass                   
                kwarg.update({ Attr: Dat })
            elif obj is Oarray:
                Oarray_fields.update({Attr: [obj, dim, init, dat_type]})
        newD = DataObj(self.nN, self.nT, [Oarray_fields], **kwarg)
        newD.parse_data( verbose )
        M =  Model( DataObj = newD ) 
        M._make_special_selections()  
        newD.MODEL[:] = M
        return M     


    def newModel_Cat(self, selObj_List, axis='T', inc_keys=[], verbose=True):
        """Takes itself and one or a list of objects and concatenates them either along the N axis (adding atoms to a frame)
        or along the T axis, adding time frames to a model"""        
        from DataObj import DataObj, Oarray, Darray
                
        #create a list of self and all objects to add
        if type(selObj_List) in SelObj_types:
            selObj_List = [self, selObj_List]
        elif type(selObj_List) in [Oarray, list, array]:
            selObj_List = [self] + list(selObj_List)            
            
        #check if dimensions are compatible (if dim=T --> all should have the same number of Atoms (N) and vice versa)
        if axis == 'T':
            if not array([ (o.nN == self.nN) for o in selObj_List ]).all():
                raise AssertionError( "The objects have different number of Atoms (%s)! and can therefore not be concatenated in time"%"".join([ str(o.nN)+', ' for o in selObj_List]) )
            T_total = array([o.nT for o in selObj_List]).sum()
            N_total = self.nN            
                  
        elif axis == 'N':
            if not array([ (o.nT == self.nT) for o in selObj_List ]).all():
                raise AssertionError( "The objects have different number of time frames (%s)! and can therefore not be concatenated"%"".join([ str(o.nT)+', ' for o in selObj_List]) )
            T_total = self.nT
            N_total = array([o.nN for o in selObj_List]).sum() 
                    
        else:
            raise ValueError("valid dimensions are only 'T' or 'N'!")

        # the actual concatenating machine
        kwarg = {}
        Oarray_fields = {}
        
        for Attr, [obj, dim, init, dat_type] in self.D.attributes.items():
            Attr = Attr.upper()

            # for all Data Objects concatinate the data (and increment the numbers if nescessary) and put them in the dict kwarg
            if obj is Darray:
                try:
                    ax = dim.index(axis) # see if we have to concatenate along this axis
                    D_concat = []        # a list of all Data arrays to be concatenated
                    for selObj in selObj_List:
                        if not hasattr( selObj.D, Attr ): 
                            raise AssertionError("All objects must have the same Data properties!")
                        
                        d = getattr( selObj, Attr.lower()+'_uf' ) # use the unflattend version 
                        if Attr in inc_keys and (len(D_concat) != 0):
                            d += D_concat[-1].max() + 1
                            
                        D_concat.append( d ) # get a copy of the Data array and put it into the list                  
                    
                    #TODO find a better way to get the Data there!                       
                    DD = concatenate( tuple(D_concat), axis=ax ) # will give an automatic copy
                    D = Darray( DD.shape, dim, DD.dtype)
                    D[:] = DD                     
                    
                except ValueError: # if it does not have the axis to concatenate along just take the data from the self obj
                    D = getattr( selObj, Attr.lower()+"_uf" ).copy() # make a copy of it to be independent
                
                kwarg.update({ Attr : D })
                
            # for all object properties we have to recreate them so we just store their existence
            elif obj is Oarray:
                if not Oarray_fields.has_key(Attr):
                    Oarray_fields.update({Attr: [obj, dim, init, dat_type]})   
        
        newD = DataObj(N_total, T_total, [Oarray_fields], **kwarg)
        newD.parse_data( verbose )#verbose=false!
        M =  Model( DataObj = newD ) 
        M._make_special_selections()  
        newD.MODEL[:] = M
        return M        
 
    def SliceT(self, T_index, verbose=True):
        """Retruns a Model of the same Selection but with a new DataObj ('.D') behind it containing only the atoms in the selection.
            Use this to reduce the number of atoms prior to reading in a trajectory!"""
        try:
            nnT = len(arange(self.nT)[T_index])
            if nnT == 0 : raise Exception                
        except:
            raise AssertionError("T_index must sclice the timeframes --> True/False array of length T or an array/list of valid indices")
        
        from DataObj import DataObj, Oarray, Darray        
        kwarg = {} #get all attributes of the old D object
        Oarray_fields = {}
        for Attr, [obj, dim, init, dat_type] in self.D.attributes.items():
            Attr = Attr.upper()
            if obj is Darray:
                Dat = getattr( self.D, Attr )
                Dat = Dat[ self.__slice_unflat(Dat.dshape) ]
                try:
                    T_dim   = Dat.dshape.index('T')
                    T_slice = tuple( [slice(None)]*T_dim + [T_index])
                    Dat = Dat[T_slice]
                except:
                    pass                
                kwarg.update({ Attr: Dat })
            elif obj is Oarray:
                Oarray_fields.update({Attr: [obj, dim, init, dat_type]})
        newD = DataObj(self.nN, nnT, [Oarray_fields], **kwarg)
        newD.parse_data( verbose )
        M =  Model( DataObj = newD ) 
        M._make_special_selections()  
        newD.MODEL[:] = M
        return M


    def combineModel(self, selObj_List, inc=True, verbose=True):
        """takes itself plus a selection Objects (or a list of selectionobjects into a new (Data) Model and returns it)
        inc=True --> will increment the atom numbers and residue numbers with each added object"""
        from DataObj import DataObj, Oarray, Darray
        
        #create a list of self and all objects to add
        if type(selObj_List) in SelObj_types:
            selObj_List = [self, selObj_List]
        elif type(selObj_List) in [Oarray, list, array]:
            selObj_List = [self] + list(selObj_List)
 
        # compatibiliy checks:
        for selObj in selObj_List:
            if(self.nT != selObj.nT): raise AssertionError("The two objects must have the same number of Time frames!")

        #new Atom size
        N_total = array([selObj.nN for selObj in selObj_List]).sum()
        
 
        kwarg = {}
        Oarray_fields = {}
        
        # these keys will be incremented to make unique atoms
        inc_keys = ['I_NDX', 'I_MODEL', 'I_RESIDUE', 'I_ATOM']
        
        for Attr, [obj, dim, init, dat_type] in self.D.attributes.items():
            Attr = Attr.upper()

            # for all Data Objects concatinate the data (and increment the numbers if nescessary) and put them in the dict kwarg
            if obj is Darray:                
                i_dat = []
                for selObj in  selObj_List:
                    if not hasattr( selObj.D, Attr ): raise AssertionError("All objects must have the same Data properties!")
                    
                    #see if we want to increment the atom number etc.
                    i = getattr( selObj.D, Attr )[:]                                                                                           
                    if(Attr in inc_keys) and (len(i_dat) != 0): i += i_dat[-1].max() - i.min() + 1 #make the numbering fit for the next part
                    i_dat.append(i)
                
                #we have to copy over from an array because r_[ returns an array not a Darray... #TODO find a better way of doing this 
                a = concatenate(i_dat)
                D = Darray(a.shape, getattr(self.D , Attr).dshape, dtype=a.dtype, buffer=a )
                #D = a[:]               
                kwarg.update({ Attr : D })               
                            
                
            # for all object properties we have to recreate them so we just store their existence
            elif obj is Oarray:
                if not hasattr( selObj.D, Attr ): raise AssertionError("All objects must have the same Object properties!")
                if not Oarray_fields.has_key(Attr):
                    Oarray_fields.update({Attr: [obj, dim, init, dat_type]})        
                    
       
        newD = DataObj(N_total, self.nT, [Oarray_fields], **kwarg)
        newD.parse_data( verbose )#verbose=false!
        M =  Model( DataObj = newD ) 
        M._make_special_selections()  
        newD.MODEL[:] = M
        return M
             

#### Coordinates geometry operation
    def cartesian_outlines(self):
        """ returns the lowest xyz coordiantes and the highest xyz coordinates of the given selection """
        return (self.x.min( self.xNaxis ),  self.x.max( self.xNaxis ))    

    def inside_box(self, box_min, box_max):
        """ Returns all atoms of a selection that are inside the cartesian rectangle
            spanned by the two pints box_min and box_max """
        b_min = (self.x > box_min).all(1)
        b_max = (self.x < box_max).all(1)
        mm = b_min * b_max
        return SelObj(  self.D, self.N.crt[mm], self.T.crt  )

    def overlapping(self, selObj, mindist=0.5):
        """Returns a selection with the Atoms that are found to be overlapping
           below a distance of mindist (in [nm]) two atoms of the respective selections are considered overlapping  
           TODO: see how we do it with multiple timeframes...."""
        if self.unity == 'A': mindist *= 10
        box_min, box_max    = selObj.cartesian_outlines() 
        R                   = self.inside_cartesian_box(box_min, box_max)
        # squared distance matrix
        OL = (R.x - selObj.x[:,newaxis])**2
        OL = OL.sum(2).min(0) < (mindist**2)        
        return R.N[ OL ]
    
    def inside_cartesian_box(self, box_min, box_max):
        """ Returns all atoms of a selection that are inside the cartesian rectangle
            spanned by the two pints box_min and box_max """
        b_min = (self.x > box_min).all(1)
        b_max = (self.x < box_max).all(1)
        mm = b_min * b_max
        return SelObj(  self.D, self.N.crt[mm], self.T.crt  )

    def inside_scew_box(self, center, box_vector_matrix):
        pass

    def inside_cyliner(self, center, axis, radius, radius_neg=None):
        pass

    def inside_sphere(self, center, radius):
        pass
  
    def cog(self):
        "center of geometry"
        return self.x.mean(self.xNaxis)     

    def com(self, weights=None):
        """center of mass,
        * weight (optional): provide own weights that will be used instead of the masses"""
        
        if self.nN == 1:    #for single Atom
            return self.x        
        
        if weights is None:
            if not hasattr(self, "m"):
                raise UserWarning("no masses found using equal weights")
                m = ones(len(self))
            else:
                m = self.m
                if (self.m == 1.).all():
                    raise UserWarning( "all Masses are one: check if you have read in correct masses!" )                
        else: # if weights are specified
            try:
                if len(weights) != len(self): raise Exception
            except:
                AssertionError("If you want to give weights provide as many as there are Atoms! (Nb atoms:%d ;Weights given:%d)"%(self.nN,len(weights)))
            m = weights
            
        # sum( Coord * weight ) / sum(weight) 
        return ( self.x * m[:,newaxis] ).sum( self.xNaxis ) / m.sum()


    def com_residue(self, weights=None):
        if self.nN == 1:    #for single Atom
            return self.x        
        return array([ r.T[self.T.crt].com() for r in self.residue ])
    
    def com_pbc(self, weights=None):
        "get the center of mass and take into account the periodic bondaries (only for rectangular boxes...)"
        if self.nN == 1:    #for single Atom
            return self.x

        if weights is None:
            if not hasattr(self, "m"):
                raise UserWarning("no masses found using equal weights")
                m = ones(len(self))
            else:
                m = self.m
                if (self.m == 1.).all():
                    raise UserWarning( "all Masses are one: check if you have read in correct masses!" )                
        else: # if weights are specified
            try:
                if len(weights) != len(self): raise Exception
            except:
                AssertionError("If you want to give weights provide as many as there are Atoms! (Nb atoms:%d ;Weights given:%d)"%(self.nN,len(weights)))
            m = weights
            
        ### the actual calculation!       
        Box_full = self.box.sum( self.axN ) # vector from the origin to the diagonal box edge
        Box_half = 0.5*Box_full                # vector to the half box
        C_Normal = self.x                                                       # the non moved coordinates
        C_Moved  = (self.x + Box_half[ self.onN ]) % Box_full[ self.onN ]         # the coordinates moved by half a box        
        S_STD    = C_Normal.std( self.axN ) < C_Moved.std( self.axN )     # (Tx3) vector that indicates in which regime the std was lower        
        M_Normal = ( C_Normal * m[ self.onN, newaxis ] ).sum(self.axN) / m.sum()               # calculating the com in both regimes 
        M_Moved  = ( C_Moved  * m[ self.onN, newaxis ] ).sum(self.axN) / m.sum() + Box_half        
        return where(S_STD, M_Normal, M_Moved)                                # select the result with the lowest std


    def _Quaternion_largest_EigenVector(self, Reference, Fitgroup=None):
        """Hidden Function used in RMSD calculation as well as fitting
            calculates the largest eigenvector between the fitgroup and the reference"""
        from numpy.linalg import eig 
        R = Reference;
        F = Fitgroup
        
        if R.nT == 1 and F.nT == 1: #only one structure in both cases                
            IP = ( (R.x - R.com()).T * (F.x - F.com()).T[:,newaxis] ).sum(2)
                
        elif R.nT == 1 and F.nT != 1: # one template multiple fitting structures
            #try:
            IP = ( (R.x - R.com()).T[newaxis,:,:,newaxis] * (F.x - F.com()[:,newaxis]).swapaxes(0,2)[:,newaxis] ).sum(2)
            #except MemoryError: #slower but less memory intensive
            #    RC = R.x - R.com()
            #    IP = array([ (RC.T * F.T[i].x.T[:,newaxis]).sum(2) for i in range(F.nT) ]) # TODO
                
        elif R.nT !=1 and  F.nT != 1: # multiple templates, multiple fitting structures
            #try:
            IP = ( (R.x - R.com()[:,newaxis]).swapaxes(0,2) * (F.x - F.com()[:,newaxis]).swapaxes(0,2)[:,newaxis] ).sum(2)
            #except MemoryError: #slower but less memory intensive
            #    IP = array([ ( (R.T[i].x - R.T[i].com()).T  * F.T[i].x.T[:,newaxis]).sum(2) for i in range(F.nT) ])
                            
        ## construct K Matrices for each time (4x4xT)
        x=0; y=1; z=2
        K = empty((4,4,F.nT)) # make space
        #diagonal elements
        K[0,0] = IP[x,x] + IP[y,y] + IP[z,z]
        K[1,1] = IP[x,x] - IP[y,y] - IP[z,z]
        K[2,2] =-IP[x,x] + IP[y,y] - IP[z,z]
        K[3,3] =-IP[x,x] - IP[y,y] + IP[z,z]
        #off diagonal (symmetric)
        K[1,0] = K[0,1] = IP[z,y]-IP[y,z]
        K[2,0] = K[0,2] = IP[x,z]-IP[z,x]
        K[3,0] = K[0,3] = IP[y,x]-IP[x,y]
        K[1,2] = K[2,1] = IP[y,x]+IP[x,y]
        K[1,3] = K[3,1] = IP[z,x]+IP[x,z]
        K[2,3] = K[3,2] = IP[z,y]+IP[y,z]            
        ## get the eigenvector with the highest eigenvalue
        RVec = empty((4,F.nT))
        RVal = empty((F.nT))
        for i in range(F.nT):
            EVal, EVec = eig( K[:,:,i] )
            RVec[:,i] = EVec[ : , argsort( EVal )[-1] ]
            RVal[i]   = EVal[ argsort( EVal )[-1] ]            
        return RVec, RVal # largest eigenvector of quarternion



    def RMSD(self, Reference):
        """ Calculates the RMSD value between two structures 
            Reference    --> a selection object with the same number of Atoms and one frame or as many frames as the current selection            
        """
        R = Reference;
        S = self;
        
        if not(R.nN == S.nN):
            raise AssertionError("Reference must have the same number of Atoms as the current selection (Reference: %d, Selection: %d)"%(R.nN, S.nN))        
        if not R.nT in [ 1, S.nT ]:
            raise AssertionError("reference must have either one frame or as many as the Selection (Reference: %d, Selection: %d)"%(R.nT, S.nT))       

        RVec, RVal = self._Quaternion_largest_EigenVector(R, S) # get the largest Eigenvector from Quarternion matrix
        
        #inner Product or Reference
        if R.nT == 1:
            Gr =  (( R.x - R.com() )**2).sum()
        else:
            Gr =  (( R.x - R.com()[:,newaxis])**2).sum(1).sum(1) #sum over all but the Time axis
            
        #inner Product or the Selection
        if S.nT == 1:
            Gs =  (( S.x - S.com() )**2).sum()
        else:
            Gs =  (( S.x - S.com()[:,newaxis])**2).sum(1).sum(1) #sum over all but the Time axis
            
        S = (Gr + Gs - 2*RVal) / S.nN
        S[ S<0 ] = 0            
        RMSD = sqrt(S)
        return RMSD
    
    def fit(self, Reference, Fitgroup=None, rot=True, trans=True, weight="mass"):
        """ Fits itself to one reference structure or a series of reference structures of the same number of frames.
            Reference    --> a selection object with the same number of Atoms and one frame or as many frames as The current selection
            Fitgroup     --> the selection of which the fit will be used to change the coordinates of the selection (if None (default) it will use the selection itself)
            rot          --> rotate the two objects to align them
            trans        --> translate the two objects to align them
            weight       --> vector of the length N (#of atoms) to weigh in the fit
            Method taken "from http://journals.iucr.org/a/issues/2005/04/00/sh5029/index.html" Theobald et al.
        """       
        from numpy.linalg import eig 
        
        #sanity checks
        if not rot and not trans:
            raise AssertionError("You have to either fit translational (trans=True) or rotational (rot=True) or both otherwise fitting is futile")       
        
        R = Reference;
        #handel default Values 
            # Fitgroup
        if Fitgroup is None: F = self
        else:                F = Fitgroup
            #weight
        if weight in ['mass','m','Mass','MASS','M']: w = R.m
        elif weight is None:                         w = ones(R.nN)
        else:
            if not weight.shape is R.m.shape:
                raise AssertionError("the weights must be a vector with the same length as the number of Atoms (%d)"%R.nN)       
            else:
                w = weight
                
        ###Check if the Fitgroup has as many frames as the object itself 
        if not(self.nT == F.nT):
            raise AssertionError("Fitgroup must have the same number of time frames as the current Object (Fitgroup: %d, current Obj: %d)"%(R.nT, self.nT))
        ##Check if two objects are compatible
        if not(R.nN == F.nN):
            raise AssertionError("Reference must have the same number of Atoms as the Fitgroup (Reference: %d, Fitgroup: %d)"%(R.nN, F.nN))        
        if not R.nT in [ 1, F.nT ]:
            raise AssertionError("reference must have either one frame or as many as the Fitgroup (Reference: %d, Fitgroup: %d)"%(R.nT, F.nT))       
       
        #store old Center of masses before Fits
        COM_R = R.com()
        COM_F = F.com()
        #center around zero        
        self.translate( -F.com() )
        if rot: # see if we have to calculate the rotation matrix
            RV, RVal = self._Quaternion_largest_EigenVector(R, F) # get the largest Eigenvector from Quarternion matrix                         
            RV2 = RV**2
            #construct rotation Matrix
            w=0; x=1; y=2; z=3 #quarternion ....
            RMat = array([ [ 1 - 2*(RV2[y]  + RV2[z]),       2*(RV[x]*RV[y] - RV[z]*RV[w]),    2*(RV[x]*RV[z] + RV[y]*RV[w])  ],
                           [ 2*(RV[x]*RV[y] + RV[z]*RV[w]),  1 - 2*(RV2[x]  + RV2[z]),         2*(RV[y]*RV[z] - RV[x]*RV[w])  ],
                           [ 2*(RV[x]*RV[z] - RV[y]*RV[w]),  2*(RV[y]*RV[z] + RV[x]*RV[w]),    1 - 2*(RV2[x]  + RV2[y])       ]])
            ## do the rotation!
            for i, SO in enumerate(self.T):
                SO.x = dot( SO.x, RMat[:,:,i] )         
        ## if translational fit was used...
        if trans: #... move to center of mass of fitgroup
            self.translate( COM_R )       
        else:# ... else move back to own center of mass
            self.translate( COM_F )                    
        return self
     
    def main_axis(self, EigenVal=False, verbose=False):
        """Returns the main axis of the Selection as three Vectors sorted ascending (Tx3x3)
           This may be useful to find the main axis of an alpha helix or a full protein.
           it is acheved by generating the covariances of all atoms for a time frame and calculating the eigenvectors"""
        from numpy.linalg import eig
        if self.nN == 1: raise ValueError("can't compute main axis from a single Atom!")
        MainAxis = empty( (self.nT,3,3), dtype = float32) #reserve Memory
        if EigenVal: EV = empty( (self.nT,3), dtype = float32)
        if verbose: print("Calculating main axis for %d Atoms in %d Time frames..."%(self.nN, self.nT))
        for i, SO in enumerate(self.T):                             # loop over timeframes!
            if verbose: print('\r...frame %d of %d'%(i+1, self.nT))
            evalue, evec = eig(cov(SO.x, rowvar=0))
            asort = argsort(evalue)[::-1]
            MainAxis[i,] = evec.T[asort]
            if EigenVal: EV[i,] = evalue[asort]
        if EigenVal:
            return MainAxis, EV  
        else:
            return MainAxis

    def pbc(self, start=0, stop=None, skip=1, verbose=True):
        """Makes all Coordinates of the selection fit into the given box for all time frames indicated"""
        from numpy.linalg import inv
        if verbose: print "Setting all Atoms into the box, ",        
        #we go to the D object to avoid cases of one (3x3) or mutltiple (Tx3x3) time frames and check if all off diagonal elements are zero
        if ( (self.D.BOX[start:stop:skip,0,1:3]  ==0).all() and  
             (self.D.BOX[start:stop:skip,1,0:3:2]==0).all() and
             (self.D.BOX[start:stop:skip,2,0:2]  ==0).all() ):
            
            #We have a rectangular box --> faster algorythm (p % diag(B))
            if verbose: " assuming rectangular box..."        
            for i, SO in enumerate(self.T[start:stop:skip].T):                
                if verbose: print "\r\t... frame %i"%i,                
                SO.x = SO.x % diag(SO.box)         
            
        else: # we have to do some matrix multiplications B * ((p(N) * B^-1) % 1)
            if verbose: "using algorithm for scew boxes..."                         
            for i, SO in enumerate(self.T[start:stop:skip].T):            
                if verbose: print "\r\t... frame %i"%i,
                SO.x = dot(   dot( SO.x, inv(SO.box) ) % 1,  SO.box   ) 
                      
        if verbose: print "...done"
        return self
        
    def pbc_mol(self, start=0, stop=None, skip=1, verbose=True, cut=2.):
        """ correct the pbc to make residues whole over the boundary (only for square Boxes so far...)"""  
        self.pbc()
        if self.unity == 'nm': cut /= 10.0
        Box_diag_vect    = self.box.sum(0)
        # get all atoms that are two bond length of the lower surfaces of the Box
        Border     = (-self.inside_box( ones(3)*cut,  Box_diag_vect - ones(3)*cut  )).by_residue()
        for dim in range(3): #for X,Y and Z
            #get the full residues that are involved
            Resi = Border.N[ Border.x[:,dim] < cut ].by_residue()
            #every atom whose coordinates are greater that half the box length at the current dimension
            High = Resi.N[ Resi.x[:,dim] > Box_diag_vect[dim]/2.0  ]
            Resi = High.by_residue()
            Low  = Resi.N[ Resi.x[:,dim]  < Box_diag_vect[dim]/2.0 ]
            # decides which part of the cut residue to move to the other side 
            CX = [ r.com()[dim] for r in Low.residue] < Box_diag_vect[dim]/2.0
            # now move the part of the reside that is to be moved  
            (Resi.residue[-CX].asSelection() & Low ).translate(   self.box[dim])
            (Resi.residue[ CX].asSelection() & High).translate(  -self.box[dim])
            
        return self
    
    def pbc_center(self, coords, start=0, stop=None, skip=1, verbose=True, pbc_mol=False, center_timefix = False):
        """Centers the box onto the given coordinates  
           coords           --> may be one coordinate for all frames (3) or a coordinate for all frames (Tx3)
           pbc_mol          --> does the pbc correction with whole residues
           center_timefix   --> takes the average box center for all time frames
        
        """
        t = len(self.time[start:stop:skip])
        if not( coords.shape == (3,) or coords.shape == (t,3) ):
            raise AttributeError("the coordinates supplied must be either a single Vector (3) or a vector for each time frame (T[start:stop:skip]x3) Vector")        
        box_center = self.box[start:stop:skip].sum(1)/2.0   #coordinates of the box center
        
        if center_timefix: # average the box center over all frames
            box_center = box_center.mean(0) 
                  
        self.translate( box_center - coords )
        if pbc_mol:
            self.pbc_mol(start, stop, skip, verbose)
        else:            
            self.pbc(start, stop, skip, verbose)

    def translate(self, vector, axis='T'):
        "translate all atoms by vector"
        if vector.shape == (3,):
            self.x += vector        
        elif axis == 'T' and len(vector) == self.nT:
            self.x += vector[:,newaxis]
        elif axis == 'T' and len(vector) == self.nN:
            if self.T == 1:
                self.x += vector[newaxis] 
            else:
                self.x += vector         

    def rotate(self, R_matrix=None, Axis=None, Axis2=None, Angle=None, Angle_in_rad=False, orig=None, verbose=False, return_R_matrix=False):
        """Rotation of all atoms in the selection:
           supply either:
           -> a rotation matrix (3x3) (R_matrix)
           -> an axis to rotate around by an angle (Axis and Angle)
           -> two axis where rotation will be performed like rotating one axis onto the other (Axis and Axis2)
           orig specifies the point to rotate around and will be the center of mass if not specified
           return_R_matrix = True will return the rotation matrix used for rotation instead of actually rotating the selection"""
        
        def _R_from2Axis(A1, A2):
            A1 = A1.T; A1 /= sqrt( (A1**2).sum(0) ); A1 = A1.T # normalize
            A2 = A2.T; A2 /= sqrt( (A2**2).sum(0) ); A2 = A2.T # normalize
            Angle = arccos( (A1*A2).sum(1) )                   # getting the angle between the vectors
            Axis  = cross( A1, A2 )                            # getting the axis perpendicular to the two vectors 
            return _R_Matrix( Axis, Angle)                     # feed it into the routine to get the Rotation Matrix
        
        def _R_fromAxisAngle(A1, Angle, A1_is_normed=False):
            if not A1_is_normed:
                A1 = A1.T; A1 /= sqrt( (A1**2).sum(0) ); A1 = A1.T # normalize 
            return _R_Matrix( Axis, Angle)
           
        def _R_Matrix(Ax_normed, Angle_rad):
            cp = cos(Angle_rad); a_sp = sin(Angle_rad)*Ax_normed
            R_matrix = array([  [ cp,         -a_sp[2],    a_sp[1] ],
                                [ a_sp[2],     cp,        -a_sp[0] ],
                                [-a_sp[1],     a_sp[0],    cp      ]]) + (1-cp) * outer(Ax_normed,Ax_normed)
            return R_matrix
            
        def _Rotate(R_matrix):
            self.translate(-orig) #move Origin to (0,0,0)
            if verbose: print("Rotating selection of %d atoms and %d time frames..."%(self.nN, self.nT))
            if R_matrix.shape == (self.nT,3,3):                         
                for i, R, SO in enumerate( zip(R_matrix, self.T) ):
                    if verbose: print '\r...frame %d of %d'%(i+1, self.nT),   
                    SO.x = dot(SO.x, R)
            elif R_matrix.shape == (3,3):
                R = R_matrix
                for i, SO in enumerate( self.T ):
                    if verbose: print '\r...frame %d of %d'%(i+1, self.nT),   
                    SO.x = dot(SO.x, R)
            else:
                raise AttributeError("Rotation matrix %s does not fit the data %s"%(R_matrix.shape, self.x.shape))
            print()
            self.translate( orig) #move back
                    
        if orig is None:
            orig = self.com() #use the center of mass (of every time frame) as the rotation origin

        if (R_matrix is not None) and (Axis is None) and (Axis2 is None) and (Angle is None): #Roatation Matix
            try:
                if R_matrix.shape != (3,3) and R_matrix.shape != (self.T,3,3): raise TypeError()
            except:
                raise AttributeError("R_matix must be a 3x3 Array or a Tx3x3 Array")

        elif (R_matrix is None) and (Axis is not None) and (Axis2 is None) and (Angle is not None): # Rotate around axis
            if not( shape(Axis) == (3,) ) or ( shape(Axis) == (self.nT,3) ):
                raise AttributeError('Axis 1 has to be either a vector (3) or a vector for every time frame (Tx3)! here it is %s'%Axis.shape )
            if not Angle_in_rad:
                Angle = deg2rad(Angle)
            R_matrix = _R_fromAxisAngle(Axis, Angle)
              
        elif (R_matrix is None) and (Axis is not None) and (Axis2 is None) and (Angle is not None): # rotate like rotating Axis onto Axis2 the other
            if not( shape(Axis) == (3,) or shape(Axis) == (self.nT,3) ):
                raise AttributeError('Axis 1 has to be either a vector (3) or a vector for every time frame (Tx3)! here it is %s'%Axis.shape )
            if not( shape(Axis2) == (3,) or shape(Axis2) == (self.nT,3) ):
                raise AttributeError('Axis 2 has to be either a vector (3) or a vector for every time frame (Tx3)! here it is %s'%Axis2.shape )
            
        else:
            raise AttributeError("""To rotate, you must either give:
            1) a rotation matrix (either 3X3 or Tx3x3)     [R_Matrix = ... ]
            2) a an axis (3 or Tx3), and angle (1 or T)    [Axis= ... , Angle= ...]
            3) two axis (each 3 or Tx3)                    [Axis= ... , Axis2= ... ]""")         
            
        if return_R_matrix:
            return R_matrix
        else:
            _Rotate(R_matrix)           


    def perm_reduction(self, Reference=None, pbc=True):
        """will relabel the residues to best fit the reference frame"""
        from NGMX_geometry import linAssignment
             
        #
        if Reference == None:       # if no reference is given use the first frame
            Reference = self.T[0]
            
        if not self.has_same_NAtoms(Reference) :            raise AssertionError("Reference must have the same number of atoms")
        if not self.residue.size == Reference.residue.size: raise AssertionError("Reference must have the same number of residues")
        if not Reference.residue.all_same_NAtoms() :        raise AssertionError("Reference must consist of only same residue type")
        if not self.residue.all_same_NAtoms() :             raise AssertionError("Selection must consist of only same residue type")
                
        if Reference.nT != 1: raise AssertionError("reference can only be one frame")       
        Rcom = Reference.com_residue()
        Scom = self.com_residue()       
        N    = Reference.residue[0].nN
        

        for SO, scom in zip(self.T, Scom):
            D = scom - Rcom[:,newaxis]
            
            if pbc: #so far only square boxes are supported                
                T = D > SO.box.diagonal()/2 * SO.box.diagonal()/2
                D = D-T
                
            D = (D**2).sum(2)
            residue_order = linAssignment(D)        
        
            atom_order = concatenate([ ((residue_order*N)+n)[:,newaxis] for n in range(N) ], axis=1).flatten()
            SO.x = SO.x[atom_order]
            

        
        
        

    #### writing out Data
    def write(self, Structure_File, start=0, stop=-1, skip=1, Title="numpymacs Structure"):
        "write the selection to a file"
        from File_IO import KnownFiles
        if type(Structure_File) is str:
            try:
                ext = os.path.splitext(Structure_File)[-1]
                Structure_File = KnownFiles[ext](Structure_File)
            except:
                raise IOError("File %s can't be used as a structure file"%Structure_File)            
        elif type(Structure_File) not in KnownFiles.values():
            raise IOError("File %s can't be used as a structure file"%Structure_File)    
        Structure_File.write_selection(self, start, stop, skip, Title)

    
class Model(SelObj): ###################################################################################################################################
    def __init__(self, structure_filename=None, trajectory_filename=None, index='all', start=1, stop=-1, skip=1, DataObj=None, verbose=True):
        if DataObj is None:
            from DataObj import DataObj
            D = DataObj(1,1) # empty dataobject
        else:            
            D = DataObj
        SelObj.__init__(self, D, arange(D.N))        
        self.load(structure_filename, trajectory_filename, index, start, stop, skip) #load structure and Trajectory (if given ...)            
            
    def load(self, structure_filename=None, trajectory_filename=None, index='all', start=1, stop=-1, skip=1):
        if structure_filename is not None:
            if not os.path.isfile(structure_filename):
                raise IOError("File %s not found"%structure_filename)
            self.load_Structure(structure_filename, index)
            self.Ncrt = arange(self.D.N)
            self.Dcrt = slice(None)
            self.D.MODEL[:] = self        
            SelObj.__init__(self, self.D, self.Ncrt, self.Dcrt)            
            self._make_special_selections()
        if trajectory_filename is not None:
            if not os.path.isfile(trajectory_filename):
                raise IOError("File %s not found"%trajectory_filename)
            self.load_Trajectory(trajectory_filename, start, stop, skip)
        return self
    
    def load_Structure(self, Structure_File, index='all'):
        from File_IO import KnownFiles
        from DataObj import DataObj
        if type(Structure_File) is str:
            try:
                ext = os.path.splitext(Structure_File)[-1]
                self.Structure_File = KnownFiles[ext](Structure_File)
            except:
                raise IOError("File %s can't be used as a structure file"%Structure_File.filename)            
        elif type(Structure_File) in KnownFiles.values():
            self.Structure_File = Structure_File        
        else:
            raise IOError("File %s can't be used as a structure file"%Structure_File.filename)
        
        self.D = DataObj(self.Structure_File.N, self.Structure_File.T)
        self.Structure_File.read_structure(self.D, index=index)
        self.D.parse_data()
        return self

    def load_Trajectory(self, Trajectory_File, start=0, stop=-1, skip=1):
        #TODO chack if file exists...
        from File_IO import KnownFiles       
        if type(Trajectory_File) is str:
            try:
                ext = os.path.splitext(Trajectory_File)[-1]
                self.Trajectory_File = KnownFiles[ext](Trajectory_File)
            except:
                raise IOError("File %s can't be used as a trajectory file"%Trajectory_File.filename)
            
        elif type(Trajectory_File) in KnownFiles.values():
            self.Trajectory_File = Trajectory_File
        
        else:
                raise IOError("File %s can't be used as a trajectory file"%Trajectory_File.filename)        
        self.Trajectory_File.read_trajectory(self.D, self.D.I_NDX, start, stop, skip)
        return self

    def __repr__(self):
        return SelObj.__repr__(self, "Model")

class Chain(SelObj):###################################################################################################################################
    def __repr__(self):
        return SelObj.__repr__(self, "Chain")

class Residue(SelObj):###################################################################################################################################
    def __repr__(self):            
        return SelObj.__repr__(self, "Residue")    

class Prot_Residue(Residue):
#    phi : C(i-1) N(i) CA(i) C(i)
#   psi : N(i) CA(i) C(i) N(i+1)
#   omega: CA(i) C(i) N(i+1) CA(i+1)
    def __repr__(self):
        return SelObj.__repr__(self, "Protein Residue") 

class Atom(SelObj):###################################################################################################################################
    def __init__(self, D, crt):
        #if len(crt) != 1:
        #    raise AssertionError("Criterion for Atom initialization has more than one index...") 
        SelObj.__init__(self, D, crt)
    
    def _crit_mask(self): # overrite from Selobj
        return self.crt

    def __repr__(self):
        return "<ATOM: %s of %s: %d>"%(self.atomname, self.resname, self.i_residue)
  
SelObj_types = [SelObj, Model, Chain, Residue, Atom]  
    
    
