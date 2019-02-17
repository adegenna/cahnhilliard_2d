import numpy as np
import sys

class CahnHilliardState():
    """
    Class for tracking the state of 1D Cahn Hilliard equations.

    **Inputs**

    ----------
    C0 : numpy float array
        Initial state.
    """
    
    def __init__(self,C0):
        M,N     = C0.shape
        self.x  = np.linspace(0,1,M)
        self.y  = np.linspace(0,1,N)
        xx,yy = np.meshgrid(self.x, self.y)
        self.M  = M
        self.N  = N
        self.xx = xx
        self.yy = yy
        self.C  = self.state2D_to_1D(C0)

    def state2D_to_1D(self,C):
        return C.ravel().reshape((1,-1))

    def state1D_to_2D(self,C):
        return C.reshape([self.M,self.N])

    def push_state_to_stack(self,Cnew):
        """
        Method to push new state to stack.
        """
        Cnew        = self.state2D_to_1D(Cnew)
        self.C      = np.vstack( [self.C, Cnew] )

    def write(self,outfile,tskip=1):
        C = self.C[::tskip] # tskip period sampling
        np.savetxt(outfile,C)
