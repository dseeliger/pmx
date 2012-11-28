import os
from ctypes import cdll, util, c_char, c_char_p, POINTER, c_int, c_long, c_int16, c_uint8, c_uint, c_float, c_void_p, byref, Structure
from numpy import array, empty, int32, float32, ones
from numpy.core import asarray
from numpy.ctypeslib import ndarray, ndpointer



#Python high level Functionsdef 

#Clib_path = "/home/dkoepfe/local/lib/python/site-packages/numpymacs/C_libs/"
Clib_path = './C_libs/'

## Libraries
lib_c   = cdll.LoadLibrary("libc.so.6")

lib_gmx = util.find_library('gmx')
if (lib_gmx == None): lib_gmx = Clib_path+"gmx.so"
libgmx=cdll.LoadLibrary(lib_gmx)

lib_md  = util.find_library('md')
if (lib_md  == None): lib_md  = Clib_path + "md.so"
libmd=cdll.LoadLibrary(lib_md)

lib_xdr  = util.find_library('xdrfile')
if (lib_xdr  == None): lib_xdr  = Clib_path + "xdrfile.so"
lib_xdr=cdll.LoadLibrary(lib_xdr)

##C Datatypes
rvec_p    = ndpointer(c_float)
Box_p     = ndpointer(c_float)
c_int_p   = POINTER(c_int)
c_uint8_p = POINTER(c_uint8)
c_float_p = POINTER(c_float)
rvec=c_float*3
matrix=c_float*3*3

##C structures
class XDRFILE(Structure):
    _fields_ = [("File",     c_void_p),
                ("XDR",      c_void_p),
                ("mode",     c_char_p),
                ("buf1",     c_int_p),
                ("buf1size", c_int),
                ("buf2",     c_int_p),
                ("buf2size", c_int)]
    

###C Functions:

# get the current position of the file pointer
def get_c_file_pointer_position(fileptr):
    c_ftell = lib_c.ftell
    c_ftell.restype = c_int
    c_ftell.argtypes = [c_void_p]
    if type(fileptr) == POINTER( XDRFILE ): # if it is an XDR file pointer use the File object of that pointer
        return c_ftell(fileptr.contents.File)
    else:
        return c_ftell(fileptr)

def set_c_file_pointer_position(fileptr, byte_offset, reference_point="start"):
    c_fseek = lib_c.fseek
    c_fseek.restype = c_int
    c_fseek.argtypes = [c_void_p, c_long, c_int]
    byte_offset = int(byte_offset)    
    reference_point = {0:"start", 1:"current_pos", 2:"end"}[reference_point]
        
    if type(fileptr) == POINTER( XDRFILE ): # if it is an XDR file pointer use the File object of that pointer
        return c_fseek(fileptr.contents.File, byte_offset, reference_point)
    else:
        return c_fseek(fileptr, byte_offset, reference_point)
 
def xtc_open(filename,mode):
    "open an xdr type file (likle edr and xtc)"
    c_xdrfile_open = lib_xdr.xdrfile_open
    c_xdrfile_open.restype  = POINTER(XDRFILE)
    c_xdrfile_open.argtypes = [c_char_p, c_char_p]
    return c_xdrfile_open( filename, mode )
    
def xtc_close(fileptr):
    lib_xdr.xdrfile_close(fileptr)

def xtc_NT(fileptr):
    N = c_int()
    T = c_int()
    lib_xdr.xdrfile_NT(fileptr, byref(N), byref(T) )
    return N.value, T.value

def xtc_read_int(fileptr, N=1):
    "int xdrfile_read_int(int *ptr, int ndata, XDRFILE *xfp);"
    ptr   = (c_int*N)()
    ndata = c_int(N)    
    lib_xdr.xdrfile_read_int( byref(ptr) , ndata,  fileptr)
    #    raise IOError('could not read %s integers'%N)
    return asarray(ptr, dtype=int32)

def xtc_read_float(fileptr, N=1):
    "int xdrfile_read_int(int *ptr, int ndata, XDRFILE *xfp);"
    ptr   = (c_float*N)()
    ndata = c_int(N)    
    lib_xdr.xdrfile_read_float( byref(ptr) , ndata,  fileptr)
    #    raise IOError('could not read %s floats'%N)
    return asarray(ptr, dtype=float32)

def xtc_read_coords2(fileptr, N, T, start, stop, skip, ndx):
    "int xdrfile_decompress_coord_float(float *ptr,  int *ncoord, float *precision, XDRFILE *xfp);"
    x 		= empty( (T,N,3), dtype=float32 )
    box		= empty( (T,3,3), dtype=float32 )
    time	= empty(  T		, dtype=float32 )
    cN		= c_int(N)
    cstart  = c_int(start)
    cstop   = c_int(stop)
    cskip   = c_int(skip)    
    c_xdrfile_read_trajectory = lib_xdr.xdrfile_read_trajectory
    c_xdrfile_read_trajectory.argtypes = [ c_void_p, 						#filepointer
										   c_int,							#N
										   ndpointer(dtype=float32, flags='C_CONTIGUOUS'),        #box 
                                           ndpointer(dtype=float32, flags='C_CONTIGUOUS'),        #time
                                           ndpointer(dtype=float32, flags='C_CONTIGUOUS'),        #x
										   c_int, c_int, c_int,				                      #start, stop, skip
										   ndpointer(dtype=float32, flags='C_CONTIGUOUS')] 		  #ndx
    realT = c_xdrfile_read_trajectory(fileptr, cN, box, time, x, cstart, cstop, cskip, ndx)
    if realT < 0:
        raise IOError("XTC reading error")    
    return realT, x, box, time

#drfile_read_trajectory(XDRFILE *xfp, int N, float *box, float *time, float *x, int start, int stop, int skip, int *ndx)
#xdrfile_read_trajectory(XDRFILE *xfp, static int N, static int T, float *box, float *time, float *x, int start, int stop, int skip, int *ndx, int len_ndx)



def xtc_read_coords(fileptr, N=1):
    "int xdrfile_decompress_coord_float(float *ptr,  int *ncoord, float *precision, XDRFILE *xfp);"
    ptr       = empty( (N,3), dtype=float32  )
    ncoord    = c_int(N)
    precision = c_float()
    c_xdrfile_decompress_coord_float = lib_xdr.xdrfile_decompress_coord_float
    c_xdrfile_decompress_coord_float.argtypes = [ ndpointer(dtype=float32, flags='C_CONTIGUOUS'),
                                               c_int_p,
                                               c_float_p,
                                               c_void_p ]
    c_xdrfile_decompress_coord_float( ptr, ncoord, precision, fileptr )
    return ptr

def xtc_read_headder(fileptr):
    magic, N, step = xtc_read_int(fileptr, 3)
    if magic != 1995:
        raise IOError("Magic number error in XTC file!")
    f = xtc_read_float(fileptr, 10)
    time = f[0]; box = f[1:].reshape([3,3])    
    #X = xtc_read_coords(fileptr, N)
    return {'magic': magic, 'N':N, 'step':step, 'time':time, 'box':box}


def rmsd(atomvect1, atomsvec2, masses=None):
    if len(atomvect1) is not len(atomsvec2):
        raise IOError( "You can only compare two equally large selections here selection one has %d Atoms and two %d Atoms"%( len(atomvect1),len(atomsvec2) ) )
    N   = len(atomvect1)
    Atoms1 = atomvect1.ctypes.as_data(c_float_p)
    Atoms2 = atomsvec2.ctypes.as_data(c_float_p)
    if masses is None:
        masses = ones(N, dtype=float32) # if no masses are given weigh all atoms equally
    masses = masses.ctypes.as_data(c_float_p)
    return libgmx.calc_similar_ind(c_int(0), c_int(N), c_void_p(), masses, Atoms1, Atoms2).value


def calc_fit_R(Atoms1, Atoms2, weights=None):
    """AtomsX are Nx3 array with the coordinates
	function will calculate rotation matrix that will best rotate Atoms2 onto Atoms1 A1 = R*A2"""
    if len(Atoms1) != len(Atoms2):
        raise IOError( "You can only compare two equally large selections here selection one has %d Atoms and two %d Atoms"%( len(Atoms1),len(Atoms2) ) )
    N   = len(Atoms1)	
    R_matrix = empty((3,3), dtype=float32)
    if weights is None:
        weights = ones(N, dtype=float32) # if no masses are given weigh all atoms equally		
        c_calc_fit_R = libgmx.calc_fit_R
        c_calc_fit_R.argtypes = [ c_int, c_int,
							  ndpointer(dtype=float32, flags='C_CONTIGUOUS'), #weights
							  ndpointer(dtype=float32, flags='C_CONTIGUOUS'), #atoms1
							  ndpointer(dtype=float32, flags='C_CONTIGUOUS'), #atoms2
							  ndpointer(dtype=float32, flags='C_CONTIGUOUS'), #R_matrix
							  ]	
    libgmx.calc_fit_R(c_int(3), c_int(N), weights, Atoms1, Atoms2, R_matrix)
    return R_matrix



  

#class C_FILE(Structure):
#    "http://en.allexperts.com/q/C-1587/2008/5/FILE-Structure.htm aus stdio.h"
#    _fields_ = [("level",    c_int),      # fill/empty level of buffer
#                ("flags",    c_uint),     # File status flags
#                ("fd",       c_char),     # File descriptor
#                ("hold",     c_uint8),    # Ungetc char if no buffer
#                ("bsize",    c_int),      # Buffer size
#                ("buffer",   c_uint8_p),  # Data transfer buffer
#                ("curp",     c_int),      # Current active pointer
#                ("istemp",   c_uint),     # Temporary file indicator
#                ("token",    c_int16)]    #  Used fdef xtc_get_number_of_frames():or validity checking


















