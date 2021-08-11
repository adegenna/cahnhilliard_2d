import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class SimData:

    def __init__( self , nx , ny , nfiles , statefile='../cpp/swig/C_' ):

        self.xx , self.yy = self.make_grid( nx , ny )
        self.nfiles = nfiles
        self.solnbase = statefile
    
    @property
    def nx( self ):
        return self.xx.shape[0]

    @property
    def ny( self ):
        return self.xx.shape[1]

    def make_grid( self , nx , ny ):
        
        x     = np.arange(nx)
        y     = np.arange(ny)
        xx,yy = np.meshgrid(x,y)
        xx = xx.T; yy = yy.T

        return xx , yy

    def read_swig_soln_file( self , i  ):
        
        filenametotal = self.solnbase + str(i) + '.out'
        print( filenametotal )
        return np.genfromtxt( filenametotal ).reshape( [self.nx,self.ny] , order = 'C' )



def plot_state_snapshot( n , simdata , ax=None ):

    if ax is None:
        fig   = plt.figure( 10 , figsize=(8,8) )
        ax    = fig.gca()

    return ax.contourf( simdata.xx , simdata.yy , simdata.read_swig_soln_file( n ) , \
                        30 , vmin=-0.4 , vmax=0.4 , cmap='afmhot' )

def make_simdata_video( simdata , mp4name='simvideo.mp4' , fig=None ):

    if fig is None:
        fig   = plt.figure( 10 , figsize=(8,8) )
    ax    = fig.gca()

    frames = [ ( plot_state_snapshot( i , simdata , ax ) ).collections for i in range(int(simdata.nfiles)) ]
    
    ani = animation.ArtistAnimation( fig , frames , interval=20 , blit=True )
    ani.save('ch2d.mp4' )