import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct
import sys

class CahnHilliardSolver():
    """
    Class for computing physics of 1D Cahn Hilliard equations.

    **Inputs**

    ----------
    inputs : InputFile
        InputFile specifying various things needed for constructor.
    state : CahnHilliardState
        Initial state needed for constructor.
    """
    def __init__(self,inputs,state,physics):
        self.inputs     = inputs
        self.outdir     = inputs.outdir
        self.saveperiod = inputs.saveperiod
        self.state      = state
        self.t_steps    = inputs.t_steps
        self.physics    = physics

    def solve(self):
        """
        Method to perform time integration of CH physics.
        """
        self.state.write(self.outdir, 0)
        for i in range(1,self.t_steps+1):
            C = self.physics.compute_update()
            self.physics.update_state(C)
            if ((i % self.saveperiod) == 0):
                self.state.write(self.outdir, i)
