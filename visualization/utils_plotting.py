from typing import Any, Callable
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from utils_postproc import *
from typing import List


def plot_state_snapshot( n : int , simdata : SimData , ax=None , strTitle='' ):

    """
    strTitle : optional string for title
    """

    if ax is None:
        fig   = plt.figure( 10 , figsize=(8,8) )
        ax    = fig.gca()

    u = simdata.read_swig_soln_file( n )

    uvar = np.var( u )

    ax.set_title( simdata.work_dir + ' , ' + strTitle )

    ax.set_xticks([])
    ax.set_yticks([])

    return ax.contourf( simdata.xx , simdata.yy , u , \
                        30 , vmin=-0.4 , vmax=0.4 , cmap='afmhot' )

def make_simdata_video( simdata , mp4name='simvideo.mp4' , fig=None ):

    if fig is None:
        fig   = plt.figure( 10 , figsize=(8,8) )
    ax    = fig.gca()

    frames = [ ( plot_state_snapshot( i , simdata , ax ) ).collections for i in range(int(simdata.nfiles)) ]
    
    ani = animation.ArtistAnimation( fig , frames , interval=20 , blit=True )
    ani.save('ch2d.mp4' )


def plot_mc_results( popt : PlottingOptions ,
                     simdata_base : SimData , 
                     fig=None , 
                     f_pass : Callable = lambda x : True , 
                     f_sort_2d : Callable[ [ List[SimData] , int , int ] , List[SimData] ] = None ):
    
    """
    plots n_files chosen from n_mc total possibilities

    inputs :
        popt : instance of PlottingOptions
        simdata_base : instance of SimData
        fig : handle to matplotlib figure for plotting
        f_pass : function for downselecting MC instances
                evaluates to False if you should skip it
        f_sort_2d : function to sort 2d array of plots in x and y directions
                inputs are List[SimData] and nrows , ncols of plot figure
    """

    if fig is None:
        fig   = plt.figure( 10 , figsize=(8,8) )
    ax    = fig.gca()

    f_mc_dir = lambda idx : simdata_base.base_dir + "/" + popt.mc_workdir_name + str(idx) + "/" + simdata_base.base_filename

    simdata_list = ( SimData( simdata_base.nx , simdata_base.ny , simdata_base.nfiles , f_mc_dir(i) ) 
                                for i in np.random.choice( np.arange(popt.ntotal)+1 , popt.ntotal , replace=False ) )
    
    nrows = int( np.floor( np.sqrt(popt.nplot) ) )
    ncols = int( np.ceil( popt.nplot / nrows ) )

    idxPlot = 1
    nplots  = 0
    idxCount = 0

    simlist = []

    while ( ( nplots < popt.nplot ) & ( idxCount < popt.ntotal ) ):
    
        si = next( simdata_list )
        idxCount += 1

        u = si.read_swig_soln_file( popt.snapshot_num )
        
        if ( u is not None ) & ( f_pass( u ) ):
                    
            simlist.append( si )

            idxPlot += 1
            nplots += 1

    if f_sort_2d is not None:

        simlist = f_sort_2d( simlist , nrows , ncols )

    P1 = []
    P2 = []

    for (i,si) in enumerate( simlist ):

        p1 = si.read_param( 'params' , -1 )
        p2 = si.read_param( 'params' , 0 )

        P1.append(p1)
        P2.append(p2)

        axi = plt.subplot( nrows , ncols , i+1 )
        plot_state_snapshot( popt.snapshot_num , si , axi , '(%.2f,%.2f)' %( p1,p2 ) )

    plt.figure()
    plt.scatter( P1 , P2 , s=25 )
    plt.xlabel( r'm' , fontsize=25 )
    plt.ylabel( r'$\chi$' , fontsize=25 )
    
    plt.show()

    return ax

