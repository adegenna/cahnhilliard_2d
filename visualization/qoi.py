import numpy as np
import matplotlib.pyplot as plt

from utils_postproc import *


def compute_surface_area_of_gradient_field( u : np.ndarray ) -> np.ndarray:
    
    """
    computes |\grad(u)|^2

    inputs: 
        u : np.array of shape ( nx , ny )

    outputs:
        gradu2 : np.array of shape ( nx , ny )
    """

    gradu2 = np.zeros_like( u )
    gradu2[1:-1,1:-1] = ( u[2:,1:-1] - u[0:-2,1:-1] )**2 + ( u[1:-1,2:] - u[1:-1,0:-2] )**2

    return gradu2
    


