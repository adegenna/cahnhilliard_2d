import numpy as np
import matplotlib.pyplot as plt
import sys,os

def main():

    # Remove any temporary files for communication with the solver
    os.system( 'rm m_*.dat complete_*.dat' )

    # Set temporal profile for parameter values
    T_amp    = np.linspace(1.0,1.0,5)
    T_x      = np.linspace(0,64,5)
    T_y      = np.linspace(32,32,5)
    T_sigma  = np.linspace(32,32,5)

    dt = 0.02
    tf = 0.1
    
    # Run solver
    count = 0
    os.system( './run_jfnk.sh &' )

    # Check filesystem for indication from solver that it is waiting for next m value
    while True:
        timestamp  = '_{:0.4f}.dat'.format( (count+1) * dt )
        petsc_done = os.path.exists( 'complete' + timestamp )
        sim_done   = os.path.exists( 'complete_sim.dat' )
        if sim_done:
            print("DRIVER PROGRAM DONE")
            break
        else:
            if petsc_done:
                outlist = [ T_amp[count+1] , T_x[count+1] , T_y[count+1] , T_sigma[count+1] ]
                np.savetxt( 'T' + timestamp , outlist , fmt='%.4f' )
                count += 1
            
if __name__ == '__main__':
    main()
