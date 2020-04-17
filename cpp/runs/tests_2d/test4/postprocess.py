import numpy as np
import os,sys,subprocess

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
            elif k == '-dt_output':
                settings.dt_out = float(v.split('\n')[0])
            elif k == '-initial_temperature_file':
                settings.T0_filename = v.split('\n')[0]
            elif k == '-da_grid_x':
                settings.nx = int(v.split('\n')[0])
            elif k == '-da_grid_y':
                settings.ny = int(v.split('\n')[0])
            
    return settings

def main( petscfile ):
    
    settings  = parse_inputs_from_petscfile( petscfile )
    
    num_bin_files = 0
    for file in os.listdir("./"):
        if ( file.endswith(".bin") & file.startswith('c_') ):
            num_bin_files += 1
    
    for i in range( num_bin_files ):
        timestamp  = '_{:0.4f}.bin'.format( i * settings.dt_out )
        solnfile_i = 'c' + timestamp
        print( '../../../build/postprocess ' + petscfile + ' ' + solnfile_i )
        os.system( '../../../build/postprocess ' + petscfile + ' ' + solnfile_i )

    return



if __name__ == '__main__':
    
    main( sys.argv[1] )
