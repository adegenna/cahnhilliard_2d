import argparse
import numpy as np
import matplotlib.pyplot as plt
from cahnhilliard_2d.src.InputFile import *
from cahnhilliard_2d.src.CahnHilliardState import *
from cahnhilliard_2d.src.CahnHilliardSpectral import *
from cahnhilliard_2d.src.CahnHilliardSolver import *

import time

def main():
    """
    Main driver script for integrating the 2D Cahn Hilliard equations.

    **Inputs**

    ----------
    args : command line arguments
        Command line arguments used in shell call for this main driver script. args must have a inputfilename member that specifies the desired inputfile name.

    **Outputs**

    -------
    inputs.outdir/results.txt : text file 
        time-integrated state
    """

    # Read inputs
    parser  = argparse.ArgumentParser(description='Input filename');
    parser.add_argument('inputfilename',\
                        metavar='inputfilename',type=str,\
                        help='Filename of the input file')
    args   = parser.parse_args()
    inputs = InputFile(args);
    inputs.printInputs();

    # Read state and grid files
    C0   = np.genfromtxt(inputs.initialstatepath , delimiter=',')

    # Problem setup
    state      = CahnHilliardState(C0)
    physics    = CahnHilliardSpectral(inputs, state)
    solver     = CahnHilliardSolver(inputs,state,physics)

    # Solve
    start = time.time()
    solver.solve()
    end   = time.time()
    print("Time elapsed: %.2f seconds" %(end-start))
    
if __name__ == '__main__':
    main()
