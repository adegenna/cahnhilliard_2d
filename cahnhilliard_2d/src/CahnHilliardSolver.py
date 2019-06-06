import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct
from scipy.integrate import solve_ivp
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
        tf = (self.t_steps+1)*self.dt
        t  = np.arange(0 , tf , self.dt*self.saveperiod)
        C0 = self.state.C.ravel()
        C  = solve_ivp( self.physics.compute_RHS , [0 , tf] , C0 , t_eval = t)
        for i in range(C.y.shape[1]):
            np.savetxt(self.outdir + "C_" + str(i*self.saveperiod) + ".out",C.y[:,i] )
