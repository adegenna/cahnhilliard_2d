import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct
import sys
from abc import ABC, abstractmethod

class CahnHilliardPhysics(ABC):
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
        self.saveperiod = inputs.saveperiod
        self.state      = state
        self.xx         = state.xx
        self.yy         = state.yy
        self.t_step     = 0
        self.t_steps    = inputs.t_steps
        self.dt         = inputs.dt
        self.epsilon    = inputs.epsilon
        super().__init__()
        
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
    
    def reset_state(self):
        """
        Method used to reset state history to just the initial condition.
        """
        self.state.reset()

    @abstractmethod
    def compute_update(self):
        pass
