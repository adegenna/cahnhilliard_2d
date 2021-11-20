from utils_plotting import *
from dataclasses import dataclass
from typing import List
import glob


def f_sort_2d( simlist : List[ SimData ] , nx : int , ny : int ) -> List[ SimData ] : 

    assert( len(simlist) <= nx * ny )

    X = []
    m = []

    for si in simlist:

        pi = np.genfromtxt( glob.glob( si.base_dir + '/params_*' )[0] , delimiter=',' )[:,-1] # hacky way to read list of parameter values 
        X.append( pi[0] )
        m.append( pi[-1] )

    sortedx = np.argsort( X )

    sorted = []

    for i in range(nx): # sort rows by y
        
        idxI = sortedx[ i*ny : (i+1)*ny ]
        idxIsort = [ idxI[j] for j in np.argsort( [ m[jj] for jj in idxI ] ) ]
        sorted.append( [ simlist[j] for j in idxIsort ] )

        #sorted.append( simlist[ idxI[ np.argsort( m[idxI] ) ] ] )

    return [ item for sublist in sorted for item in sublist ]
    


def f_sort_1d( simlist : List[ SimData ] , nx : int , ny : int ) -> List[ SimData ] : 

    assert( len(simlist) <= nx * ny )

    X = []
    m = []

    for si in simlist:

        pi = np.genfromtxt( glob.glob( si.base_dir + '/params_*' )[0] , delimiter=',' )[:,-1] # hacky way to read list of parameter values 
        X.append( pi[0] )
    
    return [ simlist[j] for j in np.argsort(X) ]




def main( popt : PlottingOptions ):

    simdata = SimData( 150 , 150 , 1 , "../data/mc_M0p2/C_" )

    if popt.plot_type == 'video':
        make_simdata_video( simdata )

    elif popt.plot_type == 'snapshot':
        plot_state_snapshot( popt.snapshot_num , simdata )
        plt.show()

    elif popt.plot_type == 'mc_results':
        plot_mc_results( popt , simdata  , f_pass=lambda x : ( np.var(x) > 0.08 ) , f_sort_2d=f_sort_1d )



if __name__ == "__main__":
    
    main( PlottingOptions( nplot=16 , ntotal=128 , snapshot_num=100 ) )