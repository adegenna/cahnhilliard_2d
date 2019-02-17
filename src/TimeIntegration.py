import numpy as np
import sys

class TimeIntegration():
    """
    Class for integrating 1D Cahn Hilliard equations.

    **Inputs**

    ----------
    state : CahnHilliardState
        Initial state needed for constructor.
    physics : CahnHilliardPhysics
        Physics object needed for constructor.
    """
    
    def __init__(self,state,physics):
        self.state      = state
        self.physics    = physics

    def euler(self):
        """
        Explicit euler method.
        """
        t_steps = self.physics.t_steps
        dt      = self.physics.dt
        C       = self.state.C
        grid    = self.state.grid
        # Integrate
        for i in range(t_steps):
            t     = self.physics.t_step
            C     = self.physics.apply_BCs(C)
            self.physics.set_current_state(C)
            rhs   = self.physics.compute_RHS()
            C    += rhs*dt
            self.physics.update_state(C)
            
