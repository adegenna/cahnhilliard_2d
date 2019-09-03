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
    
    settings.T0_filename = 'initial_temperature.dat'
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
            elif k == '-da_grid_x':
                settings.nx = int(v.split('\n')[0])
            elif k == '-da_grid_y':
                settings.ny = int(v.split('\n')[0])
            
    return settings

def write_temperature_laser_parameters_for_petsc_solver( T_amp , T_x , T_y , T_sigma , filename ):

    outlist = [ T_amp , T_x , T_y , T_sigma ]
    np.savetxt( filename , outlist , fmt='%.4f' )

    return

def main():

    # Remove any temporary files for communication with the solver
    os.system( 'rm c_*.out complete_*.out' )

    settings = parse_inputs_from_petscfile( 'petscrc.dat' )

    # Set temporal profile for parameter values
    num_changes = int( settings.tf / settings.dt )
    T_amp       = np.ones( num_changes )
    th          = np.linspace( 0.75 * np.pi , 1.75*np.pi , num_changes )
    T_x         = settings.nx * ( 1 + 0.75 * np.cos( th ) )
    T_y         = settings.ny * ( 1 + 0.75 * np.sin( th ) )
    #T_x      = np.linspace( 0.0 * settings.nx , 1.5 * settings.nx , num_changes )
    #T_y      = settings.ny//2 * np.ones( num_changes )
    #T_y      = ( T_x - 0.5 * settings.nx )**2 / ( 0.4*0.4 * settings.nx )
    T_sigma  = np.linspace( settings.nx // 2 , settings.ny // 2 , num_changes )

    # Write initial temperature field to disk for petsc
    np.savetxt( settings.T0_filename , 1.0 * np.ones( settings.nx * settings.ny ) )
    
    # Run solver
    count = 0
    os.system( './run_jfnk.sh &' )

    # Check filesystem for indication from solver that it is waiting for next m value
    while True:
        timestamp  = '_{:0.4f}.out'.format( (count+1) * settings.dt )
        petsc_done = os.path.exists( 'complete' + timestamp )
        sim_done   = os.path.exists( 'complete_sim.out' )
        if sim_done:
            print("DRIVER PROGRAM DONE")
            break
        else:
            if petsc_done:
                write_temperature_laser_parameters_for_petsc_solver( T_amp[count] ,
                                                                     T_x[count] ,
                                                                     T_y[count] ,
                                                                     T_sigma[count],
                                                                     'T' + timestamp )
                count += 1
            
if __name__ == '__main__':
    main()
