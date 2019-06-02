import numpy as np
import sys
from cahnhilliard_2d.src.CahnHilliardPhysics import CahnHilliardPhysics
from cahnhilliard_2d.src.finite_diffs import *

class CahnHilliardCook(CahnHilliardPhysics):

    def __init__(self,inputs,state):
        super().__init__(inputs,state)

    def compute_deterministic_RHS(self,phi):
        t1 = -self.m*self.b * lapl(phi)
        t2 =  self.m*self.u * lapl(phi**3)
        t3 = -self.m*self.K * biharm(phi)
        return t1 + t2 + t3

    def compute_stochastic_RHS(self):
        pass

    def compute_boundary_conditions(self,phi):
        pass

    def compute_time_step(self,U,fU):
        return U + self.dt * fU
        
    def compute_update(self):
        U    = self.get_current_state()
        fU   = self.compute_deterministic_RHS(U) + self.compute_stochastic_RHS()
        U    = self.compute_time_step(U,fU)
        U    = self.compute_boundary_conditions(U)
        return U
