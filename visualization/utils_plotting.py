import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from utils_postproc import *



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


def plot_mc_results( idx_snapshot , mc_workdir_name , n_mc , simdata_base , fig=None ):

    if fig is None:
        fig   = plt.figure( 10 , figsize=(8,8) )
    ax    = fig.gca()

    f_mc_dir = lambda idx : simdata_base.base_dir + "/" + mc_workdir_name + str(idx) + "/" + simdata_base.base_filename

    simdata_list = ( SimData( simdata_base.nx , simdata_base.ny , simdata_base.nfiles , f_mc_dir(i) ) 
                                for i in np.random.choice( np.arange(n_mc)+1 , 25 ) )
    
    for (i,si) in enumerate(simdata_list):
        
        axi = plt.subplot( 5,5,i+1 )
        try:
            plot_state_snapshot( idx_snapshot , si , axi )
        except:
            pass
    
    plt.show()

    return ax

