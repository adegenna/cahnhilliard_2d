import numpy as np
import matplotlib.pyplot as plt
from scipy import misc
import os , shutil
from shutil import copyfile
import cahnhilliard as ch
from pce.src.Polynomial import *# Clone and install from github.com:adegenna/pce.git

# ***************************************************************
# Monte Carlo sampling program for 2D modified Cahn-Hilliard
# ON: temperature-dependent CH parameters
# ON: polymer-dependent CH parameters
# OFF: thermal dynamics
# ***************************************************************

def generate_legendre_signal( x , xmin , xmax , ymin , ymax , N , eps = 0.5 ):

    class I():
        def __init__( self , polynomial_order , prior_mu , prior_sigma ):
            self.polynomial_order = polynomial_order
            self.prior_mu         = prior_mu
            self.prior_sigma      = prior_sigma
    
    inputs      = I( N , 0 , 1 )
    l           = Legendre( inputs )
    y           = np.zeros_like( x )
    x_transform = 2 * ( x - xmin ) / ( xmax - xmin ) - 1
    for i in range(N):
        l_i    = l.evaluate_1d_polynomial( i , x_transform )
        a_i    = np.random.uniform( -1 , 1 )
        y     += a_i * l_i
    # Scale within limits
    y  =  ( np.random.uniform(1 , 1 + eps) * ( ymax - ymin ) * ( y - np.min(y) ) / ( np.max(y) - np.min(y) ) + ymin )
    y[ y < ymin ] = ymin
    y[ y > ymax ] = ymax
    return y

# ********* POLYMER PARAMETERS *********
Xmin     = 0.055
Xmax     = 0.5
N        = np.mean([ 200  , 2000])
L_repeat = (10**-9) * np.mean([ 20   , 80  ]) # meters
n_repeat = 15
L_omega  = n_repeat * L_repeat
L_kuhn   = (10**-9) * np.mean([ 0.5 , 3.0  ]) # meters
# **************************************

# *********** INPUTS ***********

# Setup simulation info object
info          = ch.SimInfo();
info.t0       = 0.0
info.nx       = 128
info.ny       = 128
info.dx       = 1./info.nx
info.dy       = 1./info.ny
info.bc       = 'neumann'
info.rhs_type = 'ch_thermal_no_diffusion'

# Set up grid for spatial-field quantities
nx                = int(info.nx)
xx,yy             = np.meshgrid( np.arange(0,1,1/info.nx), np.arange(0,1,1/info.nx) )

# Setup CH parameter info object
chparams              = ch.CHparamsVector( info.nx , info.ny );
chparams.b            = ch.DoubleVector(1.0    * np.ones(nx**2))
chparams.u            = ch.DoubleVector(1.0    * np.ones(nx**2))
chparams.m            = ch.DoubleVector(0.15   * np.ones(nx**2))
chparams.sigma_noise  = 0.0
chparams.eps2_min     = 0.0
chparams.eps2_max     = 1.0
chparams.sigma_min    = 0.0
chparams.sigma_max    = 1.0e10
chparams.T_min        = 0.1
chparams.T_max        = 1.0
chparams.T_const      = ch.DoubleVector( chparams.T_max  * np.ones(nx**2))
chparams.L_kuhn       = L_kuhn
chparams.N            = N
chparams.L_omega      = L_omega
chparams.X_min        = Xmin
chparams.X_max        = Xmax
chparams.compute_and_set_eps2_and_sigma_from_polymer_params( chparams.T_max , info )

# ******************************

# Define timescales
biharm_dt         = (info.dx**4) / np.max(chparams.eps_2)
diff_dt           = (info.dx**2) / np.max( [np.max(chparams.u) , np.max(chparams.b)] )
lin_dt            = 1.0 / np.max(chparams.sigma)

# Setup checkpointing in time
n_dt              = 1000
n_tsteps          = 100
info.t0           = 0
stiff_dt          = np.min([ biharm_dt , diff_dt , lin_dt ])
t                 = np.linspace(info.t0 , info.t0 + n_dt * stiff_dt , n_tsteps+1)
dt_check          = t[1]-t[0]

print( 'Biharmonic timescale dt_biharm = ' , biharm_dt )
print( 'Diffusion timescale dt_diff = ' , diff_dt , ' = ' , diff_dt/biharm_dt , ' dt_biharm')
print( 'Linear timescale dt_lin = ' , lin_dt , ' = ' , lin_dt/biharm_dt , ' dt_biharm')
print( 'Sampling interval = ' , dt_check / stiff_dt , ' dt_stiff' )

mc_samples = 100
for i in range(50,mc_samples):

    print( '\n************** MC SAMPLE ' + str(i+1) + ' **************\n' )
    
    # Generate random temperature signal
    outdir = './mc_' + str(i+1) + '/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    T       = generate_legendre_signal( t , t[0] , t[-1] , chparams.T_min , chparams.T_max , 10 )
    T[-10:] = np.linspace(T[-10] , chparams.T_min , 10 ) # Final quench
    np.savetxt( outdir + 'T_mc_' + str(i+1) , T )
    
    # Run solver
    info.iter        = 0
    for j in range(n_tsteps):
        info.t0          = t[j]
        info.tf          = t[j+1]
        chparams.T_const = ch.DoubleVector(T[j] * np.ones(nx**2))
        print( 't0 = ', t[j]/lin_dt, ' dt_lin , tf = ', t[j+1]/lin_dt, ' dt_lin' )
        ch.run_ch_solver(chparams,info);

    # Copy simulation data to outdir
    for basename in os.listdir('./'):
        if basename.endswith('.out'):
            pathname = os.path.join('./', basename)
            if os.path.isfile(pathname):
                shutil.move(pathname, outdir)
