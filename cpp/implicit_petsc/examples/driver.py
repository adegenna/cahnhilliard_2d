import numpy as np
import matplotlib.pyplot as plt
import sys,os

def main():

    # Remove any temporary files for communication with the solver
    os.system( 'rm m_*.dat complete_*.dat' )

    # Set temporal profile for m values
    m  = np.linspace(-0.15,0.15,10)
    dt = 0.10
    tf = 1.0
    
    # Run solver
    count = 0
    os.system( './run_jfnk.sh &' )

    # Check filesystem for indication from solver that it is waiting for next m value
    while True:
        timestamp  = '_{:0.2f}.dat'.format( (count+1) * dt )
        petsc_done = os.path.exists( 'complete' + timestamp )
        sim_done   = os.path.exists( 'complete_sim.dat' )
        if sim_done:
            print("DRIVER PROGRAM DONE")
            break
        else:
            if petsc_done:
                np.savetxt( 'm' + timestamp , m[count+1].reshape([1,1]) , fmt='%.2f' )
                count += 1
            
if __name__ == '__main__':
    main()
