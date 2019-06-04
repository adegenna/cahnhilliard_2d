import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct,diff
import sys
from cahnhilliard_2d.src.CahnHilliardPhysics import CahnHilliardPhysics

def dct2d(x,inverse=False):
    t    = 2 if not inverse else 3
    temp = dct(x,type=t,norm='ortho').transpose()
    return dct(temp,type=t,norm='ortho').transpose()

def spectral_diff2d(x):
    temp = np.array( [diff(xi) for xi in x] ).transpose()
    temp = np.array( [diff(ti) for ti in temp] ).transpose()
    return temp

class CahnHilliardSpectral(CahnHilliardPhysics):
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
        super().__init__(inputs,state)
        self.compute_laplacian_eigenvalues()
        self.compute_current_state_dct()
        
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

    def compute_stochastic_term(self,A=1):
        eta       = np.random.normal(0,1,[self.state.N , self.state.M])
        eta[0]    = 0
        eta[-1]   = 0
        eta[:,0]  = 0
        eta[:,-1] = 0
        #kx        = np.fft.fftfreq(self.state.N)*self.state.N
        #ky        = np.fft.fftfreq(self.state.M)*self.state.M
        #kkx,kky   = np.meshgrid(kx,ky)
        eta       = A * spectral_diff2d( eta )
        return eta
        
    def compute_update(self):
        """
        Method to compute an update to the CH physics.
        """
        U       = self.get_current_state()
        # compute the shifted nonlinear term
        fU      = (U*U*U) - ((1+self.inputs.eyre_a)*U)
        # compute the right hand side in tranform space
        hat_rhs  = self.hat_U + (self.Seig*dct2d(fU))
        hat_rhs += self.compute_stochastic_term(self.noise_amplitude)
        # compute the updated solution in tranform space
        self.hat_U   = hat_rhs/self.CHeig
        # invert the cosine transform
        U       = dct2d(self.hat_U,inverse=True)
        return U
