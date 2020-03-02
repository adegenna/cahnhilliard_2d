import numpy as np
#import matplotlib.pyplot as plt
import sys,os
#import h5py

class PetscSettings:

    def __init__(self):
        self.tf          = None
        self.dt          = None
        self.T0_filename = None
        self.nx          = None
        self.ny          = None
        self.nz          = None

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
            elif k == '-da_grid_z':
                settings.nz = int(v.split('\n')[0])
            
    return settings

def write_temperature_laser_parameters_for_petsc_solver( T_amp , T_x , T_y , T_sigma , filename ):

    outlist = [ T_amp , T_x , T_y , T_sigma ]
    np.savetxt( filename , outlist , fmt='%.4f' )

    return

def generate_quarter_circle_laser_path( amp , sigma_temp , nx , ny , num_changes ):

    T_amp       = amp * np.ones( num_changes )
    th          = np.linspace( 0.75 * np.pi , 1.75*np.pi , num_changes )
    T_x         = nx * ( 1 + 0.75 * np.cos( th ) )
    T_y         = ny * ( 1 + 0.75 * np.sin( th ) )
    T_sigma     = sigma_temp * np.ones( num_changes )
    
    return T_amp , T_x , T_y , T_sigma

def generate_const_global_temperature( amp , nx , ny , nz , num_changes ):

    T_amp   = amp   * np.ones( num_changes )
    T_x     = nx//2 * np.ones( num_changes )
    T_y     = ny//2 * np.ones( num_changes )
    T_z     = nz//2 * np.ones( num_changes )
    T_sigma = nx*ny*nz * np.ones( num_changes )
    
    return T_amp , T_x , T_y , T_z , T_sigma

def generate_gaussian_temperature_profile_on_3d_face( amp , sigma_temp , nx , ny , nz , x0 , y0 ):

    x        = np.arange( nx )
    y        = np.arange( ny )
    z        = np.arange( nz )
    xx,yy,zz = np.meshgrid( x , y , z )
    T_xyz    = amp * np.exp( -0.5 * ( ( xx-x0 )**2 + (yy - y0)**2 ) / sigma_temp**2 ) * np.heaviside( -(zz-1) , 1 )
    
    return xx,yy,zz,T_xyz

def generate_gaussian_temperature_profile_on_2d_face( amp , sigma_temp , nx , ny , y0 ):

    x        = np.arange( nx )
    y        = np.arange( ny )
    xx,yy    = np.meshgrid( x , y )
    T_xy     = amp * np.exp( -0.5 * ( (yy - y0)**2 ) / sigma_temp**2 ) * np.heaviside( -(xx-1) , 1 )
    
    return xx,yy,T_xy

def generate_constant_temperature_profile_2d( amp , nx , ny ):
    
    x        = np.arange( nx )
    y        = np.arange( ny )
    xx,yy    = np.meshgrid( x , y )
    T_xy     = 0. * xx + amp
    
    return xx,yy,T_xy

def main():

    builddir = '../../build/'

    settings = parse_inputs_from_petscfile( 'petscrc.dat' )

    # Set temporal profile for parameter values
    num_changes                 = int( settings.tf / settings.dt )
    sigma_temp                  = 1.0*settings.nx / 4.
    amp_temp                    = 1.0
    T_amp , T_x , T_y , T_sigma = generate_quarter_circle_laser_path( amp_temp , sigma_temp , settings.nx , settings.ny , num_changes )
    
    # Write initial temperature field to disk for petsc
    const_T             = 0.3
    initial_T           = const_T * np.ones( settings.nx * settings.ny )
    initial_T_file      = 'initial_temperature.ascii'
    np.savetxt( initial_T_file , initial_T , fmt='%1.8f' )
    os.system( builddir + 'preprocess petscrc.dat ' + initial_T_file )
    
    # Write the initial temperature source field to disk for petsc
    T_source      = 0.0 * initial_T
    T_source_file = 'initial_temperature_source.ascii'
    np.savetxt( T_source_file , T_source , fmt='%1.8f' )
    os.system( builddir + 'preprocess petscrc.dat ' + T_source_file )
    
    # Write the initial solution field to disk for petsc
    initial_U      = 0.005 * ( 2.0 * np.random.uniform(0,1,settings.nx*settings.ny) - 1.0 )
    initial_U_file = 'initial_soln.ascii'
    np.savetxt( initial_U_file , initial_U , fmt='%1.8f' )
    os.system( builddir + 'preprocess petscrc.dat ' + initial_U_file )

    # Write temperature dirichlet field to disk for petsc
    y0                  = settings.ny // 2
    amp_temp            = 1.0
    xx,yy,initial_T     = generate_gaussian_temperature_profile_on_2d_face( amp_temp , sigma_temp , settings.nx , settings.ny , y0 )
    initial_T           = initial_T.ravel()
    initial_T           = np.maximum( initial_T , 0.3 )
    initial_T_file      = 'dirichlet_T.ascii'
    np.savetxt( initial_T_file , initial_T , fmt='%1.8f' )
    os.system( builddir + 'preprocess petscrc.dat ' + initial_T_file )

    # Write ch dirichlet field to disk for petsc
    y0                  = settings.ny // 2
    amp_ch              = 1.0
    xx,yy,initial_T     = generate_constant_temperature_profile_2d( amp_ch , settings.nx , settings.ny )
    initial_T           = initial_T.ravel()
    initial_T_file      = 'dirichlet_ch.ascii'
    np.savetxt( initial_T_file , initial_T , fmt='%1.8f' )
    os.system( builddir + 'preprocess petscrc.dat ' + initial_T_file )

    # Run solver
    count = 0
    os.system( './run_test.sh &' )

    # Check filesystem for indication from solver that it is waiting for next temperature profile
    while True:
        timestamp        = '_{:0.4f}.bin'.format( (count+1) * settings.dt )
        timestamp_final  = '_{:0.4f}.bin'.format( settings.tf )
        petsc_done = os.path.exists( 'complete' + timestamp )
        sim_done   = os.path.exists( 'c' + timestamp_final )
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
