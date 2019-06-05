import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct
from scipy.integrate import odeint
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
        self.dt         = inputs.dt
        self.physics    = physics

    def solve_eyre(self):
        """
        Method to perform time integration of CH physics.
        """
        self.state.write(self.outdir, 0)
        for i in range(1,self.t_steps+1):
            C = self.physics.compute_update()
            self.physics.update_state(C)
            if ((i % self.saveperiod) == 0):
                self.state.write(self.outdir, i)

    def solve(self):
        self.state.write(self.outdir, 0)
        t  = np.arange(0 , (self.t_steps+1)*self.dt , self.dt*self.saveperiod)
        C0 = self.state.C.ravel()
        C  = odeint(self.physics.compute_RHS , C0 , t)
        print(C.shape)
        for i in range(C.shape[0]):
            np.savetxt(self.outdir + "C_" + str(i*self.saveperiod) + ".out",C[i] )
