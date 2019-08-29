import numpy as np
import matplotlib.pyplot as plt
import sys,os

class PetscSettings:

    def __init__(self):
        self.tf          = None
        self.dt          = None
        self.T0_filename = None
        self.nx          = None
        self.ny          = None

def parse_inputs_from_petscfile( petscfile ):

    settings = PetscSettings()
    
    T0_filename = 'initial_temperature.dat'
    with open( petscfile ) as f:
        petscsettings = f.readlines()
        for line in petscsettings:
            try:
                k,v = line.split(' ')
            except:
                k = None; v = None
            if   k == '-t_final':
                settings.tf = float(v.split('\n')[0])
            elif k == '-dt_check':
                settings.dt = float(v.split('\n')[0])
            elif k == '-initial_temperature_file':
                settings.T0_filename = v.split('\n')[0]
            elif k == '-nx':
                settings.nx = int(v.split('\n')[0])
            elif k == '-ny':
                settings.ny = int(v.split('\n')[0])
            
    return settings

def main():

    # Remove any temporary files for communication with the solver
    os.system( 'rm m_*.out complete_*.out' )

    # Set temporal profile for parameter values
    T_amp    = np.linspace(1.0,1.0,5)
    T_x      = np.linspace(0,64,5)
    T_y      = np.linspace(32,32,5)
    T_sigma  = np.linspace(32,32,5)

    settings = parse_inputs_from_petscfile( 'petscrc.dat' )

    # Write initial temperature field to disk for petsc
    np.savetxt( settings.T0_filename , 1.0 * np.ones( settings.nx * settings.ny ) )
    
    # Run solver
    count = 0
    os.system( './run_jfnk.sh &' )

    # Check filesystem for indication from solver that it is waiting for next m value
    while True:
        timestamp  = '_{:0.4f}.out'.format( (count+1) * dt )
        petsc_done = os.path.exists( 'complete' + timestamp )
        sim_done   = os.path.exists( 'complete_sim.out' )
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
