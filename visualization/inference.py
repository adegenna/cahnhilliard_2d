import numpy as np
import matplotlib.pyplot as plt

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import WhiteKernel, DotProduct, RBF , Matern , ConstantKernel

from qoi import *
from utils_postproc import *


def fit_gp( npts : int , 
            basedir : str = '../data/mc_M0p2' ) -> GaussianProcessRegressor:

    c = []
    p = []

    count = 0

    for i,di in enumerate( glob.glob( basedir + '/mc_*' ) ):

        s = SimData( 150 , 150 , 1 , di + '/C_' )

        c.append( s.read_swig_soln_file( 100 ) )

        p.append( s.read_params_as_dict( 'params_'  ) )

        count +=1
        if count >= npts:
            break

    gradc = [ np.sum( compute_surface_area_of_gradient_field( ci ) ) for ci in c ]

    kernel = ConstantKernel( constant_value=1e5 , constant_value_bounds=( 1e4, 1e7) ) * \
             Matern( length_scale=1e-1 , length_scale_bounds=(1e-2, 1e0) , nu=1.5 ) + \
             WhiteKernel( noise_level=1e4 , noise_level_bounds=(1e4,1e4) )
    # kernel = ConstantKernel( constant_value=1e4 , constant_value_bounds=( 1e-8, 1e8) ) * \
    #          Matern( length_scale=1e0 , length_scale_bounds=(1e-3, 1e2) , nu=1.5 ) + \
    #          WhiteKernel( noise_level=1e3 , noise_level_bounds=(1e0,1e6) )
    gpr  = GaussianProcessRegressor( kernel=kernel , \
                                     n_restarts_optimizer=100 , \
                                     #alpha=1e2 , \
                                     normalize_y=True )

    gpr.fit( np.array( [ pi['X'] for pi in p ] ).reshape(-1,1) , np.array( gradc ).reshape(-1,1) )

    return gpr


def plot_gp( gpr : GaussianProcessRegressor , 
             color : str = 'b' ):

    x_test = np.linspace( 0 , 0.7 , 1000 ).reshape(-1,1)

    y_hat, y_sigma = gpr.predict( x_test , return_std=True )

    plt.scatter( gpr.X_train_ , gpr.y_train_ + gpr._y_train_mean , c=color )

    lower = y_hat - y_sigma.reshape(-1,1)
    upper = y_hat + y_sigma.reshape(-1,1)
    plt.gca().plot( x_test , y_hat , c=color )
    plt.gca().fill_between( x_test.ravel() , lower.ravel() , upper.ravel() , alpha=0.5 )

    plt.xlabel( r'$\chi$' )
    plt.ylabel( r'$y = |\nabla u|^2 \;\; , \;\; \rho( y | \chi )$' )




def main():

    npts = [ 64 , 8 ]

    gpr = [ fit_gp(ni) for ni in npts ]

    colors = [ 'b' , 'r' ]#, 'g' ]

    for (gpi,ci) in zip( gpr , colors ):
        plot_gp( gpi , ci )

    plt.show()




if __name__ == '__main__':

    main()