# -*- coding: utf-8 -*-
"""
Modified: 29.6.15

@author: Guglielmo Saggiorato <g.saggiorato@fz-juelich.de> or <astyonax@gmail.com>
"""

import numpy as np

def pma(data,axis=0):
    """
    principal components analysis of 2D data - expects time on the 0th axis
    returns:
    sorted eigvals, and eigvectors[:,mode-id], correlation matrix
    """
    i=data[0]
    
    #MD=np.einsum('i,j->ij',i,i.T)
    MD=np.outer(i,i)
   
    for i in data[1:]:    
        MD+=np.outer(i,i)
    MD/=(data.shape[0]-1)
    
    # this may be faster to do
    # but beware of the transposition!!
    #   datat=data.T
    #   MD=np.mean(datat[:,np.newaxis]*datat,axis=-1,ddof=1)
    
    #print "computing eigenvalues"
    l,e=np.linalg.eigh(MD)
    
	# sort the eingen values and modes by ascending order of eigen value
    idxx=np.argsort(l)[::-1]
    l=l[idxx]
    e=e[:,idxx]
    return l,e,MD
    
def get_XY(data,e):
    """
    Compute time evolution of the modes amplitude 
    as the scalar product between the data and each mode 
    """
    
    #return [(data*e[:,j]).sum(axis=1) for j in xrange(e.shape[1])]

def get_XY_fast(data,e):
    """
    Compute time evolution of the modes amplitude 
    as the scalar product between the data and each mode 
    This uses a faster shortcut, but it not tested as much as the slower version
    """
    # this is not tested enough, but faster. In case of errors, switch
    # the comment btw the next two lines
    return e.T.dot(data.T)
    #return [(data*e[:,j]).sum(axis=1) for j in xrange(e.shape[1])]

def reconstruc_frommodes(e,ampl):
    ampl=np.array(ampl)
    e=np.array(e)
    
    rec=np.zeros((ampl.shape[1],e.shape[0]),dtype=np.float64)
   
    N=ampl.shape[1]
    e=e.T
    
    for eigvect,X in zip(e,ampl):
        
        rec0=np.tile(eigvect,(N,1))
        rec+=rec0*X[:,np.newaxis]
    return rec
