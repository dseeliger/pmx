from numpy import *
from numpy.linalg import norm, solve

def dihedral(A1, A2, A3, A4):#(N,3)
    "calculates the dihedral angle (for all times provided) from A1--A2-!-A3--A4 around axis between A2 and A3"
    ##Todo sanity check to hybe (Tx3) Arrays
    b1 = A1-A2; b2 = A2-A3; b3 = A3-A4
    b2_norm = sqrt((b2**2).sum(1))
    cross_12 = cross( b1, b2 )
    cross_23 = cross( b2, b3 )		
    return arctan2( b2_norm * (b1*cross_23).sum(1), (cross_12*cross_23).sum(1) ) *180 /pi

def angle(A1, A2, A3):#(N,3)
    "calculates the angle between A1--A2--A3"
    b1 = (A1-A2);   b2 = (A3-A2)
    b1 /= sqrt((b1**2).sum(1))[:,newaxis]
    b2 /= sqrt((b2**2).sum(1))[:,newaxis]
    return arctan( (b1*b2).sum(1) ) *180/pi

def dist(A1, A2):#(N,3)
    "Returns the distance between A and B for every time"
    return sqrt(((A1-A2)**2).sum(1))	

def dist_sq(A1, A2):#(N,3)
    "Returns the square of the distance (faster!)"
    return ((A1-A2)**2).sum(1)

def axis_projection(center, axis, points): #(N,3)/(3,)
    """ returns the distance of the projection from the center onto the axis as well as the distance from each point to the axis 
    Solves D = dfa*B + doa*A;
    D vector between the center and each point
    A unit vector in axis direction
    B unit difference vector between the center and each atom
    dfa = distance from axis
    doa = distance on axis 
    
    TODO --> !
    center: point in space (3,) --> (3,) (T,3)
    axis:   axis on which to project (3,) --> (3,) (T,3)
    points: number of atom coordinates (N,3) --> (3,) (T,N,3)
    
    returns array with [dfa,doa] (N,2) --> (T,N,2)
    """
    if len(points.shape) == 1:	#Compatibility for Points (N,3) or (3,)
        points = points[newaxis,:]

    #constructing vectors
    A      = axis.copy()                           # axis-vector to project onto                (3,)
    D      = points - center                       # vector connecting between center and point (N,3) 
    B      = cross(cross(A, D), A)                 # vector perprndicular to the axis pointing towards the point (N,3)
           
    #normalizing vectors so dfa and doa are in the right units (nm or A)  
    A     /= norm(A)
    B     /= sqrt((B**2).sum(1))[:,newaxis]
           
    #solve for
    i0 = A.nonzero()[0][0]       # avoid picking an axis that is zero!
    i1 = (i0+1)%3
    aa = A[i1]/A[i0]
        
    dfa = (D[:,i1]-D[:,i0]*aa) / (B[:,i1]-B[:,i0]*aa)   #distance from axis   
    doa = (D[:,i0] - dfa*B[:,i0]) / A[i0]              #distance on axis
    #return
    return c_[doa, dfa]                              #return a vector (2,N) with "distance distance on axis", "distance from axis"
     
    # calculating distances	
    #return array([ solve(c_[axis,bb,cc],pxx)[:2] for bb,cc,pxx in zip(b,c,px) ]) #if too slow consider C wrapper....
