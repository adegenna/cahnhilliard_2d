import numpy as np
import matplotlib.pyplot as plt
import sys,os
import h5py

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

def generate_quarter_circle_laser_path( amp , nx , ny , num_changes ):

    T_amp       = amp * np.ones( num_changes )
    th          = np.linspace( 0.75 * np.pi , 1.75*np.pi , num_changes )
    T_x         = nx * ( 1 + 0.75 * np.cos( th ) )
    T_y         = ny * ( 1 + 0.75 * np.sin( th ) )
    T_sigma     = np.linspace( nx // 2 , ny // 2 , num_changes )
    
    return T_amp , T_x , T_y , T_sigma

def generate_const_global_temperature( amp , nx , ny , num_changes ):

    T_amp   = amp   * np.ones( num_changes )
    T_x     = nx//2 * np.ones( num_changes )
    T_y     = ny//2 * np.ones( num_changes )
    T_sigma = nx*ny * np.ones( num_changes )
    
    return T_amp , T_x , T_y , T_sigma

def main():

    # Remove any temporary files for communication with the solver
    os.system( 'mkdir old' )
    os.system( 'mv *.out old/' )

    settings = parse_inputs_from_petscfile( 'petscrc.dat' )

    # Set temporal profile for parameter values
    num_changes                 = int( settings.tf / settings.dt )
    #T_amp , T_x , T_y , T_sigma = generate_quarter_circle_laser_path( 1.0 , settings.nx , settings.ny , num_changes )
    T_amp , T_x , T_y , T_sigma = generate_const_global_temperature( 0.5 , settings.nx , settings.ny , num_changes )

    # Write initial temperature field to disk for petsc
    #np.savetxt( settings.T0_filename , 1.0 * np.ones( settings.nx * settings.ny ) )
    with h5py.File( settings.T0_filename , 'w' ) as f:
        dset = f.create_dataset( "default" , data=1.0*np.ones( settings.nx * settings.ny ) )

    
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
