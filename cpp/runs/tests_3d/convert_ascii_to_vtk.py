import numpy as np
import matplotlib.pyplot as plt
import os,sys
from pyevtk.hl import gridToVTK

def read_mpi_soln_file( statefile , n ):
    
    c = np.zeros( n )
    count = 0
    with open( statefile ) as f:
        for line in f:
            try:
                c[count] = line
                count += 1
            except:
                pass
    
    return c

def main( statefile , nx=100 , ny=100 , nz=100 ):

    print( nx , ny , nz )
    
    c = read_mpi_soln_file( statefile , nx*ny*nz )
    c = c.reshape( (nx,ny,nz) )
    
    x = np.arange(0, nx+1)
    y = np.arange(0, ny+1)
    z = np.arange(0, nz+1)
    cname = '.'.join( statefile.split( '.' )[:-1] )

    gridToVTK( cname , x , y , z , cellData = {cname: c} )

    return

if __name__ == '__main__':
    main( sys.argv[1] ,  int(sys.argv[2]) , int(sys.argv[3]) , int(sys.argv[4]) )
