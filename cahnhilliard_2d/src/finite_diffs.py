import numpy as np

def get_internal_stencil(u):
    assert( len(u.shape) == 2)
    nx,ny = u.shape
    ul    = u[
    return ul,ur,ut,ub,uc

    
def lapl(u,dx=1.):
    assert( len(u.shape) == 2)
    nx,ny = u.shape
    uxx   = np.zeros_like(u)
    uyy   = np.zeros_like(y)
    
    

def biharm(U):
    pass
