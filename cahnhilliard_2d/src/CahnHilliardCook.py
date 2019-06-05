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
        self.dx = self.state.x[1]-self.state.x[0]
        self.dy = self.state.y[1]-self.state.y[0]

    def compute_deterministic_RHS(self,phi):
        t1 = -self.m*self.b * lapl(phi,self.dx)
        t2 =  self.m*self.u * lapl(phi**3,self.dx)
        t3 = -self.m*self.K**2 * biharm(phi,self.dx)
        return t1 + t2 + t3

    def compute_stochastic_RHS(self):
        pass

    def compute_boundary_conditions(self,phi):
        #phi = apply_neumann_bc_d1_and_d3(phi)
        phi = apply_periodic_bc(phi)
        return phi
    
    def compute_RHS(self,U,t):
        U    = self.state.state1D_to_2D(U)
        U    = self.compute_boundary_conditions(U)
        fU   = self.compute_deterministic_RHS(U) #+ self.compute_stochastic_RHS()
        fU   = self.state.state2D_to_1D(fU).flatten()
        return fU
