from numpy import *
import library as lib
import geometry as geo

class Darray(ndarray):
    "Made for Data"
    def __new__(subtype, shape, dshape='N', dtype=float, buffer=None, offset=0, strides=None, order=None):
        if type(shape) is ndarray:
            data = shape
            shape = data.shape      
            obj = ndarray.__new__(subtype, shape, dtype, buffer, offset, strides,order)
            obj[:] = data
        else:
            obj = ndarray.__new__(subtype, shape, dtype, buffer, offset, strides,order)
        if dshape in ['N','T'] or type(dshape) is int or type(dshape) is tuple:
            obj.dshape = dshape
        else: 
            raise TypeError("dshape must be a tuple, 'N', 'T' or a number!")        
        return obj

    def __array_finalize__(self,obj):
        self.dshape = getattr(obj, 'dshape', None)

class Oarray(ndarray):
    "Made for pointers to selections"
    def __new__(subtype, shape, dshape='N', dtype=object, buffer=None, offset=0, strides=None, order=None):
        obj = ndarray.__new__(subtype, shape, dtype, buffer, offset, strides,order)
        if type(dshape) in [ str, int, tuple]: # or type(dshape) is int or type(dshape) is tuple:
            obj.dshape = dshape
        else: 
            raise TypeError("dshape must be a tuple, 'N', 'T' or a number!")  
        return obj

    def __array_finalize__(self,obj):
        self.dshape = getattr(obj, 'dshape', None)
        
    ### info functions
    
    def all_same_NAtoms(self):
        """returns true if all residues have same number of atoms"""
        return len(  unique( array([r.nN for r in self]) )  ) == 1
            
        
    def asSelection(self):
        if len(self) == 1:
            return self[0]
        else:
            from SelObj import SelObj
            D_test = array([ SO.D for SO in self ])
            if not (D_test[1:] == D_test[:-1]).all(): #if not all Objects share the same Dataobject
                raise MemoryError('Not all selesctions share the same Data Object!')                
            ndx = []
            for SO in self:
                ndx += SO.N.crt.tolist()
            ndx = unique(array(ndx))         
            return SelObj( D_test[0], ndx, self[0].T.crt ) #TODO find consensus for the time frame criterion...

    def com(self):
        A = array([ s.com() for s in self ])
        if self[0].nT > 1:
            A.swapaxes(0,1)

        return Darray()

        
    def inside_cart_box(self, a, b ):
        """ returns the objects which coms are inside the cartesian box made up by the two coordinates 
            b,a = (3,) Vector with two points that span a box"""
        X = self.com()
        a,b = c_[a,b].min(1), c_[a,b].max(1)        
        return X[ (X>a) * (X<b) ]
        
    
    def inside_box(self, center, box):
        """ Retruns all elements that are (with their center of mass inside a box
            box is boxvectors (3,3) --> (x,y,z)
            center center of the box 
        """    
        from numpy.linalg import inv        
        X = self.com()
        center -= (box*0.5).sum(1) #move the center to the lower corner of the box
        P = dot( X , inv(box) )
        return self[ (P>0) * (P<1) ]
    
        
    def inside_cylinder(self, center, axis, radius, axis_length=None, ax_length_negative=None):
        #TODO
        X = self.com()        
        if axis_length == None: axis_length = (axis**2).sum()
        pass


    
    def split(self, nparts=2):
        pass
#        """ splits a selection in parts
#            if nparts is a number it will try to split in n equal parts
#            if it is an array or a list it will sort accordingly """
#            
#        if type(nparts) is int:
#            if len(self) % nparts != 0:
#                raise AssertionError("can't split %d items into %d groups"%(len(self), nparts))
#            
#        elif type(nparts) in [ndarray, list]:
#                raise AssertionError("not yet implemented... Sorry")
#            
#        else:
#            raise UserWarning("You can only give a number of equal parts or an array")
                
        
         

# attributes that are standard and are required!
attributes_req = { 'I_NDX':     [Darray,    ('N',),      arange,       uint64  ], #numbers
                   'I_MODEL':   [Darray,    ('N',),      1,            uint32  ],
                   'I_CHAIN':   [Darray,    ('N',),      1,            uint32  ],
                   'I_RESIDUE': [Darray,    ('N',),      1,            uint32  ],
                   'I_ATOM':    [Darray,    ('N',),      arange,       uint64  ],
                   
                   'ATOMNAME':  [Darray,    ('N',),      "ATOM  ",    'S6'     ], #names
                   'RESNAME':   [Darray,    ('N',),      None,        'S5'     ],
                   'CHAINNAME': [Darray,    ('N',),      None,        'S5'     ],
                                      
                   'MODEL':     [Oarray,    ('N',),      None,        object   ], #pointer arrays
                   'CHAIN':     [Oarray,    ('N',),      None,        object   ],
                   'RESIDUE':   [Oarray,    ('N',),      None,        object   ],                   
                    
                    "UNITY":    [Darray,    (1,),        'A',         'S2'     ],
                    
                    "X":        [Darray,    ('T','N',3), None ,       float32  ],
                    "BOX":      [Darray,    ('T',3,3),   1,           float32  ],
                    "TIME":     [Darray,    ('T',),      0,           float32  ],
                    "STEP":     [Darray,    ('T',),      0,           uint64   ]                   
                  }

# attributes that come into play for PDB files
attributes_pdb = {  'RACE':     [Darray,    ('N',),      "ATOM  ",    'S6'     ],
                    'BFAC':     [Darray,    ('N',),      0    ,       float32  ],
                    'OCC':      [Darray,    ('N',),      0    ,       float32  ],
                    'ALTLOC':   [Darray,    ('N',),      None ,       'S1'     ]
                  }

# attributes that help but take more time to startup and are not standard 
attributes_chem = { 'M':        [Darray,    ('N',),      1 ,          float32  ],
                    "Q":        [Darray,    ('N',),      0,           float32  ],
                    "SYMBOL":   [Darray,    ('N',),      None,        'S2'     ]
                   }  

attributes_atom = {'ATOM':      [Oarray,    ('N',),      None,        object  ],
                   }




class DataObj():
    #def __init__(self, N, T=1, attr=[attributes_req], **kwargs):
    def __init__(self, N, T=1, Sattr=[attributes_req, attributes_pdb, attributes_chem], **kwargs):
        self.T = T        # length of stored times 
        self.N = N        # length of all Data Arrays bound to atoms
        self.attributes = {}
        for key, val in kwargs.items(): 
            if type(val) in [Darray, Oarray]:
                setattr(self, key.upper() ,val)
                self.attributes.update({key.upper(): [type(val), val.dshape, None, val.dtype]})
            else:
                raise TypeError("All objects for the DataObject must either be Darrays or Oarrays! %s is of type %s"%(key, type(val))) 
               
        new_attr = {};  [ new_attr.update(d) for d in Sattr ] #other attributes are initialized here
        self.add_attributes(new_attr)

    def add_attributes(self, attr_dict):
        """Supply an attribute (An array with values) by supplying a dict with:
        name: [ObjType, dimensions, initial value, datatype]
        where
        --> name: will be converted to all uppercase for Dataobject and lowercase for the selection Objects
        --> ObjectType: "Darray" for regular data or "Oarray" for an pointer array to other Selction objects (e.g. Chains etc.)
        --> dimensions: "('N', 'T', 3)" will scale with the number of atoms, the number of time frames and have tree entries for every atom and timeframe
        --> initial value: can be a number string etc. or None to create an empty object, finally also arange is supported for onedimensional arrays
        --> datatype: can be float32 int32 or an abstract type
        example:  {"BOX": [Darray, ('T',3,3), 1, float32 ]}
        """
        for attr, [obj, dim, init, dat_type]  in attr_dict.items():
            if hasattr( self, attr.upper() ): continue # only add attributes not already given by the kwargs
            if obj in [Darray, Oarray]: # main array types Data arrays and Object arrays
                n_dim = self._attr_para_convert(attr, dim) # convert N and T to right values
                A = obj( n_dim, dim, dtype=dat_type )
                if init is arange:
                    A[:] = arange(n_dim[0])
                elif init is not None:
                    A.fill(init)
                setattr(  self, attr.upper(), A  ) 
            else:
                raise  TypeError("All objects for the DataObject must either be Darrays or Oarrays")
        self.attributes.update(attr_dict)        
            
    def _attr_para_convert(self, atr, dim):
        "converts a tuple with dimensions like ('N','T',3) --> (1000, 500, 3) for self.N=1000 and self.T=500"            
        rdim=[]        
        for d in dim:
            if d == 'N':
                rdim.append(self.N)
            elif d == 'T':
                rdim.append(self.T)
            elif type(d) is int:
                rdim.append(d)
            else:
                raise TypeError(  "Error in attribute %s: dimensionality may only contain 'N','T' or a fixed number but is %s"%( atr, str(dim) )  )
        return tuple(rdim)
        
    def set_T_size(self, newT, **kwargs):        
        if self.T == newT:
            return # nothing to do already has the right size         
        oldT = self.T
        self.T = newT
        if oldT > newT: #we want to crop the arrays 
            for attr, [obj, dim, init, dat_type]  in self.attributes.items():
                if 'T' in dim:
                    shaper = []
                    for d in dim:
                        if d == 'T': shaper.append(slice(None, newT)) # =^ XXX[:newT]
                        else: shaper.append(slice(None))
                    A = getattr( self, attr.upper() )
                    setattr( self, attr.upper(),  A[tuple(shaper)] )                             
        else: #we want to extend the arrays
            for attr, [obj, dim, init, dat_type]  in self.attributes.items(): 
                if 'T' in dim: #for all Arrays that scale somehow with T (number of timeframes)
                    shaper = []
                    for d in dim:
                        if d == 'T': shaper.append(slice(None, oldT)) # =^ XXX[:oldT]
                        else: shaper.append(slice(None))
                    n_dim = self._attr_para_convert(attr, dim)
                    if init is arange: ## init the new attribute as sepcified by the attribute
                        A = obj( arange(n_dim[0]), dim, dtype=dat_type )
                    elif init is None:
                        A = obj( n_dim, dim, dtype=dat_type )
                    else:
                        A = obj( n_dim, dim, dtype=dat_type )
                        A.fill(init)
                    A[shaper] = getattr( self, attr.upper() ) # fill array with the old data where possible
                    setattr(  self, attr.upper(), A  )
                                
    def a2nm(self):
        "convert units to Aangstroms"
        if self.UNITY == "A":
            self.X   /= 10.
            self.BOX /= 10.
            #self.V /= 10.
            self.UNITY[0] = "nm"
        elif self.UNITY[0] == "nm":
            print ("Unit is already Aangstroms")

    def nm2a(self):
        "convert units to nanometers"
        if self.UNITY[0] == "nm":
            self.X   *= 10.
            self.BOX *= 10.
            #self.V *= 10.
            self.UNITY[0] = "nm"
        elif self.UNITY[0] == "A":
            print("Unit is already nano meters")           

    def parse_data(self, verbose=True):
        "Create objects for chains, Atoms and Residues, Lookup Masses and charges ..."
        if hasattr(self, "CHAIN") or hasattr(self, "RESIDUE"):
                self._make_chains_residues(verbose)
        #if hasattr(self, "CHAIN"):
        #    self._make_chains(verbose)
        #
        #if hasattr(self, "RESIDUE"):
        #    self._make_residues(verbose)   
        
        if hasattr(self, "ATOM"):
                self._make_atoms(verbose)
                    
        if hasattr(self, "M") or hasattr(self, "SYMB"):
                self._make_symbol_masses(verbose)


    def _make_chains(self, verbose=False):
        from SelObj import Chain
        if verbose: print "\t ...chains"
        for ch in unique(self.I_CHAIN):
            CH_idx = nonzero(self.I_CHAIN == ch)[0]
            self.CHAIN[ CH_idx ] = Chain(self, CH_idx)
            
    def _make_residues(self, verbose):
        from SelObj import Residue
        if verbose: print "\t ...residues"
        for re in unique(self.I_RESIDUE):
            RE_idx = nonzero(self.RESIDUE == re)[0]
            self.CHAIN[ RE_idx ] = Residue(self, RE_idx)
    
    def _make_chainsI(self, verbose=False):
        from SelObj import Chain
        R = {}
        for i,v in enumerate(self.I_CHAIN):
            
            if R.haskey(v):     # if we already have it append index to the list
                R[v].append(i)
            else:
                R.update({ v:[i] }) # if we did not already have it, 
        
        for r in R.iterkeys():   # create the chain objects
            self.CHAIN[ r ] = Chain(self, R[r]) 

    def _make_residuesI(self, verbose):        
        from SelObj import Residue
        R = {}
        for i,v in enumerate(self.I_RESIDUE):
            if R.haskey(v):     # if we already have it append index to the list
                R[v].append(i)
            else:
                R.update({ v:[i] }) # if we did not already have it, 
        
        for r in R.iterkeys():   # create the chain objects
            self.CHAIN[ r ] = Residue(self, R[r]) 
        
    

    def _make_chains_residues(self, verbose=False):
        from SelObj import Chain, Residue 
        ## Chains
        if verbose: print "\t ...chains"
        CH_sortmask   = self.I_CHAIN.argsort()          # Create a series of indices that would sort the Chains
        CH_sorted     = self.I_CHAIN[ CH_sortmask ]     # Create a sorted array
        CH_change     = nonzero( r_[True, -( CH_sorted[1:]== CH_sorted[:-1] ), True] )[0]
        CH_indices       = [ CH_sortmask[CH_a:CH_e] for CH_a, CH_e in zip(CH_change[:-1], CH_change[1:]) ]        
        if hasattr(self, "CHAIN"):
            for CH_idx in CH_indices:
                CH_idx.sort()
                self.CHAIN[ CH_idx ] = Chain(self, CH_idx) ## here we must avoid broadcasting and store the object ...
        ## Residues 
        if hasattr(self, "RESIDUE"):
            if verbose: print "\t ...residues"    
            for CH_idx in CH_indices:                    
                RES_sortmask = self.I_RESIDUE[ CH_idx ].argsort()
                RES_sorted   = self.I_RESIDUE[ CH_idx ][ RES_sortmask ]
                RES_change   = nonzero( r_[True, -( RES_sorted[1:]== RES_sorted[:-1] ), True] )[0]
                RES_indices  = [ RES_sortmask[RES_a:RES_e] for RES_a, RES_e in zip(RES_change[:-1], RES_change[1:]) ]            
                for RES_idx in RES_indices:
                    RES_idx =  CH_idx[ RES_idx ]
                    RES_idx.sort()                
                    self.RESIDUE[ RES_idx ] =  Residue(self, RES_idx)

    def _make_atoms(self, verbose=False):
        from SelObj import Atom                
        ## Atoms
        #from Atom import Atom ## won't work without ??? don't know why
        if verbose: print "\t ...atoms"
        self.ATOM = array([ Atom(self, i) for i in range(self.N) ])

    def _make_symbol_masses(self, verbose=False):
        if verbose: print "\t...chemical symbols and atomic masses"
        Name_set = unique(self.ATOMNAME)
        for N in Name_set:            
            symb = lib.get_symbol(N)
            name_mask = where(self.ATOMNAME == N)
            if hasattr(self, "SYMBOL"):
                self.SYMBOL[name_mask]  = symb
            if hasattr(self, "M"):
                mass = 0;
                try:
                    mass = lib._atommass[symb]
                except:
                    print "no Mass found for %s. Setting it to zero"%symb
                self.M[name_mask] = mass
            
            
    
            
            
            
            
            
            
        
    