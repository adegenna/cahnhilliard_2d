import numpy as np
import sys
from cahnhilliard_2d.src.CahnHilliardPhysics import CahnHilliardPhysics
from cahnhilliard_2d.src.finite_diffs import *

class CahnHilliardCook(CahnHilliardPhysics):

    def __init__(self,inputs,state):
        super().__init__(inputs,state)
        self.m  = self.inputs.p_m
        self.b  = self.inputs.p_b
        self.u  = self.inputs.p_u
        self.K  = self.inputs.p_K
        self.sig      = self.inputs.p_sig
        self.phi_star = self.inputs.p_phi_star
        self.noise_sigma = self.inputs.p_noise_sigma
        self.boundary_conditions = self.inputs.boundary_conditions
        self.dx = self.state.x[1]-self.state.x[0]
        self.dy = self.state.y[1]-self.state.y[0]
        self.xx, self.yy = np.meshgrid(self.state.x,self.state.y)
        np.random.seed(0)
        if (self.boundary_conditions == 'neumann'):
            self.apply_boundary_conditions = apply_neumann_bc_d1_and_d3
        elif (self.boundary_conditions == 'periodic'):
            self.apply_boundary_conditions = apply_periodic_bc
        else:
            print('Error: boundary_conditions must be either neumann or periodic')
            sys.exit()

    def compute_deterministic_RHS(self,phi):
        t1 = -self.m*self.b * lapl(phi,self.dx)
        t2 =  self.m*self.u * lapl(phi**3,self.dx)
        t3 = -self.m*self.K**2 * biharm(phi,self.dx)
        t4 =  self.sig * (phi - self.phi_star)
        return t1 + t2 + t3 + t4

    def compute_stochastic_RHS(self,phi):
        nx,ny  = phi.shape
        field  = np.random.normal(0,self.noise_sigma,[nx,ny])
        field -= field.mean()
        return field
    
    def compute_stochastic_RHS_fourier(self,phi,n_wavenumbers=10):
        nx,ny = phi.shape
        field = np.zeros_like(phi)
        a     = np.random.normal(0,self.noise_sigma,[n_wavenumbers,n_wavenumbers])
        for i in range(1,n_wavenumbers):
            f_i = np.cos( i * np.pi * self.xx )
            for j in range(1,n_wavenumbers):
                f_j    = np.cos( j * np.pi * self.yy )
                field += a[i,j] * f_i * f_j
        return field

    def compute_RHS(self,t,U):
        U    = self.state.state1D_to_2D(U)
        U    = self.apply_boundary_conditions(U)
        fU   = self.compute_deterministic_RHS(U) + self.compute_stochastic_RHS(U)
        fU[0:2] = 0; fU[-2:] = 0; fU[:,0:2] = 0; fU[:,-2:] = 0;
        fU   = self.state.state2D_to_1D(fU).flatten()
        return fU
