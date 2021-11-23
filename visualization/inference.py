from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from typing import List , Tuple
from numpy.random import Generator
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import WhiteKernel, DotProduct, RBF , Matern , ConstantKernel
import pickle

from qoi import *
from utils_postproc import *

import samplers.src.methast as mcmc  # https://github.com/adegenna/samplers


class CHSolutionManager:

    def __init__( self , basedir : str = '../data/mc_M0p2' ):
    
        self._basedir = basedir
        self._c = []
        self._p = []

        for i,di in enumerate( glob.glob( basedir + '/mc_*' ) ):

            s = SimData( 150 , 150 , 1 , di + '/C_' )

            self._c.append( s.read_swig_soln_file( 100 ) )

            self._p.append( s.read_params_as_dict( 'params_'  ) )

    def make_random_subset( self , npts : int ) -> CHSolutionManager:

        """
        Make a subset manager, managing npts randomly selected data pts

        inputs : 
            npts : desired size of the random subset

        outputs : 
            CHSolutionManager to a random subset of the data managed by caller
        """
        
        subset = CHSolutionManager( self._basedir )

        idx_rand   = np.random.choice( np.arange( len(self._c) ) , npts )
        subset._c  = [ self._c[idx] for idx in idx_rand ]
        subset._p  = [ self._p[idx] for idx in idx_rand ]

        return subset

    def make_uniform_spaced_subset( self , npts : int , nameParam : str ) -> CHSolutionManager:

        """
        Make a subset manager, managing npts data points, uniformly spaced w.r.t. parameter nameParam 

        inputs : 
            npts : desired size of the random subset

        outputs : 
            CHSolutionManager to a random subset of the data managed by caller
        """

        subset = CHSolutionManager( self._basedir )

        idx_p_sort = np.argsort( [ pi[nameParam] for pi in self._p ] )[0::( int(np.floor(len(self._p)/(npts-1))) - 1 )]

        subset._c  = [ self._c[idx] for idx in idx_p_sort ]
        subset._p  = [ self._p[idx] for idx in idx_p_sort ]

        return subset

    @property
    def c( self ) -> List[ np.ndarray ]:

        """
        outputs : 
            List of solutions (each of which is a np.array)
        """

        return self._c

    @property
    def p( self ) -> List[ Dict[ str , float ] ]:

        """
        outputs : 
            List of d-dimn parameters (each is a dict mapping d string keys to d values)
        """

        return self._p



def fit_gp( npts : int , 
            basedir : str = '../data/mc_M0p2' ) -> GaussianProcessRegressor:

    cp  = CHSolutionManager( )
    sub = cp.make_uniform_spaced_subset( npts , 'X' )

    gradc = [ np.sum( compute_surface_area_of_gradient_field( ci ) ) for ci in sub.c ]

    kernel = ConstantKernel( constant_value=1e5 , constant_value_bounds=( 1e4, 1e7) ) * \
             Matern( length_scale=1e-1 , length_scale_bounds=(1e-2, 1e0) , nu=1.5 ) + \
             WhiteKernel( noise_level=1e4 , noise_level_bounds=(1e4,1e4) )
    gpr  = GaussianProcessRegressor( kernel=kernel , \
                                     n_restarts_optimizer=100 , \
                                     normalize_y=True )

    gpr.fit( np.array( [ pi['X'] for pi in sub.p ] ).reshape(-1,1) , np.array( gradc ).reshape(-1,1) )

    return gpr


def plot_gp( gpr : GaussianProcessRegressor , 
             color : str = 'b' , 
             legend : str = None ):

    x_test = np.linspace( 0 , 0.7 , 1000 ).reshape(-1,1)

    y_hat, y_sigma = gpr.predict( x_test , return_std=True )

    plt.scatter( gpr.X_train_ , gpr.y_train_ + gpr._y_train_mean , c=color )

    lower = y_hat - y_sigma.reshape(-1,1)
    upper = y_hat + y_sigma.reshape(-1,1)
    line, = plt.gca().plot( x_test , y_hat , c=color )
    line.set_label( legend )
    plt.gca().fill_between( x_test.ravel() , lower.ravel() , upper.ravel() , alpha=0.5 )
    plt.legend()

    plt.xlabel( r'$\chi$' )
    plt.ylabel( r'$y = |\nabla u|^2 \;\; , \;\; \rho( y | \chi )$' )


def plot_example_gp( ):

    npts = [ 64 , 8 ]

    gpr = [ fit_gp(ni) for ni in npts ]

    colors = [ 'b' , 'r' ]#, 'g' ]

    for (gpi,ci) in zip( gpr , colors ):
        plot_gp( gpi , ci )

    plt.show()


def load_gpr( savefile : str ) -> GaussianProcessRegressor:

    return pickle.load( open( savefile , "rb" ) )

def save_gpr( gpr : GaussianProcessRegressor , savefile : str ):

    with open( savefile , "wb" ) as f:
        pickle.dump( gpr , f )  


def main( gpr_gt_filename : str = None , gpr_coarse_filename : str = None ):

    gpr_gt = fit_gp( 128 ) if gpr_gt_filename     is None else load_gpr( gpr_gt_filename )
    gpr    = fit_gp( 4 )   if gpr_coarse_filename is None else load_gpr( gpr_coarse_filename )

    plot_gp( gpr , 'b' , 'coarse' )
    plot_gp( gpr_gt , 'r' , 'high_res' )

    X_true      = 0.4
    X_meanPrior = 0.25
    y_true , _  = gpr.predict( np.array( [ X_true ] ).reshape(-1,1) , return_std=True )
    y_scale     = 15000. / 10 # Statistical length scale for range of empirically observed y values
    X_scale     = 0.1 # Statistical length scale for range of X values

    def f_likelihood( x ):
        y_hat, y_sigma = gpr.predict( np.array( [ x ] ).reshape(-1,1) , return_std=True )
        return np.maximum( 1e-8 , np.exp( -0.5 * ( ( y_hat - y_true )**2 + (y_sigma)**2 ) / y_scale**2 ) )

    def f_prior( x ):
        return np.maximum( 1e-8 , np.exp( -0.5 * ( x - X_meanPrior )**2 / X_scale**2 ) )

    def f_target( x : List[float] ):
        return f_likelihood( x[0] ) * f_prior( x[0] )

    sampler = mcmc.Metropolis( f_target , 
                               x0={ 'X' : X_meanPrior } , 
                               k_sigma_jump=0.5 , 
                               burn_in=4096 )

    for i in range(2**15):
        _ = sampler.draw()
    
    plt.figure()

    yhist,xhist = np.histogram( [ xi['X'] for xi in sampler.samples() ] , 32 )

    plt.bar( xhist[:-1] , yhist , width=np.diff(xhist), edgecolor="black", align="edge" , label='post' )
    plt.plot( np.linspace( -1 , 1 , 1000) , f_prior( np.linspace( -1 , 1 , 1000) ) * np.max(yhist) , 'r' , label='prior' )
    plt.plot( [ X_true-1e-5 , X_true+1e-5 ] , [ 0 , 2*np.max(yhist) ] , 'k--' , linewidth=4 , label='true' )
    plt.legend()

    plt.xlim( [ 0 , 0.7 ] )
    plt.ylim( [ 0 , 1.1 * np.max(yhist) ] )

    plt.xlabel( 'X' )
    plt.ylabel( r'$\rho( X | y )$' )
    
    plt.show()


if __name__ == '__main__':

    main( 'gpr_128' , 'gpr_4' )