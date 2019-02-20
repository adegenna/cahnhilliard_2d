import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct
import sys

def dct2d(x,inverse=False):
    t    = 2 if not inverse else 3
    temp = dct(x,type=t,norm='ortho').transpose()
    return dct(temp,type=t,norm='ortho').transpose()

class CahnHilliardPhysics():
    """
    Class for computing physics of 1D Cahn Hilliard equations.

    **Inputs**

    ----------
    inputs : InputFile
        InputFile specifying various things needed for constructor.
    state : CahnHilliardState
        Initial state needed for constructor.
    """
    def __init__(self,inputs,state):
        self.inputs     = inputs
        self.outdir     = inputs.outdir
        self.saveperiod = inputs.saveperiod
        self.state      = state
        self.xx         = state.xx
        self.yy         = state.yy
        self.t_step     = 0
        self.t_steps    = inputs.t_steps
        self.dt         = inputs.dt
        self.epsilon    = inputs.epsilon
        self.compute_laplacian_eigenvalues()
        self.compute_current_state_dct()

    def set_parameters(self,p):
        """
        Method used to set parameters of physics (ie., epsilon).
        """
        self.epsilon = p
        
    def get_current_state(self):
        """
        Method to return current state (in 2D form).
        """
        C = self.state.state1D_to_2D(self.state.C)
        return C

    def update_state(self,Cnew):
        """
        Method to push new state (1D form) to stack.
        """
        Cnew         = self.state.state2D_to_1D(Cnew)
        self.state.C = Cnew

    def set_current_state(self,C):
        """
        Method to set current state (1D form).
        """
        C            = self.state.state2D_to_1D(C)
        self.state.C = C

    def compute_laplacian_eigenvalues(self):
        """
        Method to compute laplacian eigenvalue matrix and associated quantities.
        Note: assumes Neumann BCs
        """
        N  = self.state.N
        M  = self.state.M
        dx = self.state.x[1]-self.state.x[0]
        dy = self.state.y[1]-self.state.y[0]
        # time marching update parameters
        lam1 = self.inputs.dt / (dx**2)
        lam2 = self.epsilon**2 * lam1 / (dx**2)
        # unscaled eigenvalues of the laplacian (nuemann bc)
        L1   = np.tile( 2*np.cos(np.pi*np.arange(N)/(N-1)) - 2 , [M,1] ).T
        L2   = np.tile( 2*np.cos(np.pi*np.arange(M)/(M-1)) - 2 , [N,1] )
        Leig = L1 + L2
        # scaled eigenvalues of stabilized CH update matrix
        self.CHeig = np.ones((N,M)) - (self.inputs.eyre_a*lam1*Leig) + (lam2*Leig*Leig)
        # scaled eigenvalues of the laplacian
        self.Seig = lam1*Leig

    def compute_current_state_dct(self):
        """
        Method to compute dct of current state
        """
        C = self.get_current_state()
        self.hat_U = dct2d(C)

    def compute_update(self):
        """
        Method to compute an update to the CH physics.
        """
        U       = self.get_current_state()
        # compute the shifted nonlinear term
        fU      = (U*U*U) - ((1+self.inputs.eyre_a)*U)
        # compute the right hand side in tranform space
        hat_rhs = self.hat_U + (self.Seig*dct2d(fU))
        # compute the updated solution in tranform space
        self.hat_U   = hat_rhs/self.CHeig
        # invert the cosine transform
        U       = dct2d(self.hat_U,inverse=True)
        return U

    def solve(self):
        """
        Method to perform time integration of CH physics.
        """
        self.state.write(self.outdir, 0)
        for i in range(1,self.t_steps+1):
            C = self.compute_update()
            self.update_state(C)
            if ((i % self.saveperiod) == 0):
                self.state.write(self.outdir, i)
            
    def reset_state(self):
        """
        Method used to reset state history to just the initial condition.
        """
        self.state.reset()

    def reset(self):
        self.compute_laplacian_eigenvalues()
        self.compute_current_state_dct()
