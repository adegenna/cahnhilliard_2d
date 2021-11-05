import numpy as np
import matplotlib.pyplot as plt
from numpy.core.arrayprint import IntegerFormat
from scipy import misc
import cahnhilliard as ch
import csv
import os , sys
from typing import Dict, Tuple , List


def convert_temperature_to_flory_huggins( T , Tmin , Tmax , Xmin , Xmax ):
    assert( (T > 1e-5) & (Tmin > 1e-5) )
    return ch.CHparamsVector().convert_temperature_to_flory_huggins( T , Tmin , Tmax , Xmin , Xmax )

def compute_dimensionless_ch_params_from_polymer_params( X , N , L_omega , L_kuhn , m ):
    eps_2    = ch.CHparamsVector().compute_eps2_from_polymer_params( X , m , L_kuhn , L_omega , N )
    sigma    = ch.CHparamsVector().compute_sigma_from_polymer_params( X , m , L_kuhn , L_omega , N )
    return eps_2 , sigma


class MaterialParamsUniformDistribution:

    def __init__( self ):

        self.params = [ 'X' , 'N' , 'L_omega' , 'L_kuhn' , 'm' ]

        self.lower = dict( zip( self.params ,
            [ 0.055 , 
              200. , 
              ( 20e-9 ) * 15 , # [ repeat spacing in meters ] * [ number of repeats ]
              0.5e-9 , # meters
              -1.0 ]
        ) )

        self.upper = dict( zip( self.params , 
            [ 0.5 , 
              2000. , 
              ( 80e-9 ) * 15 , # [ repeat spacing in meters ] * [ number of repeats ]
              3.0e-9 , # meters
              1.0 ]
        ) )

        self.Tmin = 0.1
        self.Tmax = 1.0

    def draw( self ) -> Dict[ str , float ]:
        
        vals = np.random.uniform( list( self.lower.values() ) , list( self.upper.values() ) )

        vals[0] = convert_temperature_to_flory_huggins( 
            np.random.uniform(self.Tmin,self.Tmax) , self.Tmin , self.Tmax , self.lower['X'] , self.upper['X'] )

        return dict( zip( self.params , vals ) )

    def mean( self , strvar : str ) -> float:

        return 0.5 * ( self.lower[strvar] + self.upper[strvar] )


def setup_SimInfo( ) -> ch.SimInfo:

    info          = ch.SimInfo()

    info.t0       = 0.0
    info.nx       = 150
    info.ny       = 150
    info.dx       = 1./info.nx
    info.dy       = 1./info.ny
    info.bc       = 'neumann'

    return info


def setup_chparams( info : ch.SimInfo , p : Dict[ str , float ] ) -> ch.CHparamsVector:

    """
    inputs: 
        info : instance of ch.SimInfo
        p : dict of mc material parameters, e.g. produced by MaterialParamsUniformDistribution::draw()
    """

    eps_2 , sigma = compute_dimensionless_ch_params_from_polymer_params( 
        X=p['X'] , N=p['N'] , L_omega=p['L_omega'] , L_kuhn=p['L_kuhn'] , m=p['m'] )

    chparams      = ch.CHparamsVector( int(info.nx) , int(info.ny) )

    nx                = int(info.nx)
    xx,yy             = np.meshgrid( np.arange(0,1,1/info.nx) , np.arange(0,1,1/info.nx) )

    chparams.eps_2        = ch.DoubleVector( eps_2  * np.ones(nx**2) )
    chparams.sigma        = ch.DoubleVector( sigma  * np.ones(nx**2) )
    chparams.m            = ch.DoubleVector( p['m'] * np.ones(nx**2) )

    return chparams

def compute_linear_timescale( chparams : ch.CHparamsVector , info : ch.SimInfo ) -> float:

    biharm_dt         = (info.dx**4) / np.max(chparams.eps_2)
    diff_dt           = (info.dx**2) / np.max( [np.max(chparams.u) , np.max(chparams.b)] )
    lin_dt            = 1.0 / np.max(chparams.sigma)

    return np.min( [ biharm_dt , diff_dt , lin_dt ] )



def run_sample( chparams : ch.CHparamsVector , info : ch.SimInfo , mcparams : Dict[ str , float ] , nsample : int ):

    """
    inputs: 
        chparams : instance of ch.CHparamsVector
        info : instance of ch.SimInfo
        p : dict of mc material parameters, e.g. produced by MaterialParamsUniformDistribution::draw()
        nsample : mc sample number (e.g., this is nsample=5 of 100 total samples)
    """

    n_dt = 10

    stiff_dt = compute_linear_timescale( chparams , info )

    n_tsteps          = 20
    info.t0           = 0
    t                 = np.linspace(info.t0 , info.t0 + n_dt * stiff_dt , n_tsteps+1)

    # Run solver
    w = csv.writer( open( info.outdir + "ch_" + str(nsample) + ".csv" , "w" ) )
    w.writerow( ['Sampling interval = %.2f dt_lin = %.2e' %( ( (t[1]-t[0]) / stiff_dt) , (t[1]-t[0]) )] )

    for i in range(n_tsteps):
        info.t0          = t[i]
        info.tf          = t[i+1]
        w.writerow( ['t0 = %.2f dt_lin = %.2e , tf = %.2f dt_lin = %.2e' %( (t[i]/stiff_dt) , t[i] , (t[i+1]/stiff_dt) , t[i+1] )] )
        ch.run_ch_solver_non_thermal(chparams,info)




def main( root_outdir : str ):

    nmc  = 2
    pgen = MaterialParamsUniformDistribution()
    
    for i in range(nmc):

        outdir       = root_outdir + '/mc_' + str(i+1) + '/'
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        pi       = pgen.draw()
        info     = setup_SimInfo()
        info.outdir  = outdir
        chparams = setup_chparams( info , pi )

        run_sample( chparams , info , pi , i+1 )

        # Save mc parameters
        w = csv.writer( open( outdir + "params_" + str(i+1) + ".csv" , "w" ) )
        for key, val in pi.items():
            w.writerow( [ key , val ] )





if __name__ == '__main__':

    root_outdir = '/home/adegennaro/Projects/appmath/cahnhilliard_2d/data/mc_params'

    if len(sys.argv) != 1:
        root_outdir = sys.argv[1]
    
    main( root_outdir )