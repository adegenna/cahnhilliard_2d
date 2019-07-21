import numpy as np
import matplotlib.pyplot as plt
from scipy import misc
import cahnhilliard as ch

# ***************************************************************
# Example driver program for 2D modified Cahn-Hilliard
# ON: temperature-dependent CH parameters
# ON: polymer-dependent CH parameters
# OFF: thermal dynamics
# ***************************************************************

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
info.nx       = 100
info.ny       = 100
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
chparams.m            = ch.DoubleVector(0.0    * np.ones(nx**2))
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

print( np.sqrt(chparams.eps_2[0]) , chparams.sigma[0] )

# Define timescales
biharm_dt         = (info.dx**4) / np.max(chparams.eps_2)
diff_dt           = (info.dx**2) / np.max( [np.max(chparams.u) , np.max(chparams.b)] )
lin_dt            = 1.0 / np.max(chparams.sigma)

# Setup checkpointing in time
n_dt              = 6000
n_tsteps          = 100
info.t0           = 0
stiff_dt          = np.min([ biharm_dt , diff_dt , lin_dt ])
t                 = np.linspace(info.t0 , info.t0 + n_dt * stiff_dt , n_tsteps+1)
dt_check          = t[1]-t[0]

# Setup time-dependent temperature profile
T                            = np.zeros(n_tsteps)

# T[0:n_tsteps//4]             = chparams.T_min
# T[n_tsteps//4:n_tsteps//2]   = np.linspace( chparams.T_min , chparams.T_max , n_tsteps//4 )
# T[n_tsteps//2:3*n_tsteps//4] = chparams.T_max
# T[3*n_tsteps//4:]            = chparams.T_min

T_mid            = 0.5*( chparams.T_max + chparams.T_min )
T[0:20]          = chparams.T_max
T[20:25]         = np.linspace( chparams.T_max , T_mid , 5 )
T[25:50]         = T_mid
T[50:]           = np.linspace( T_mid , chparams.T_min , 50 )

print( 'Biharmonic timescale dt_biharm = ' , biharm_dt )
print( 'Diffusion timescale dt_diff = ' , diff_dt , ' = ' , diff_dt/biharm_dt , ' dt_biharm')
print( 'Linear timescale dt_lin = ' , lin_dt , ' = ' , lin_dt/biharm_dt , ' dt_biharm')
print( 'Sampling interval = ' , dt_check / stiff_dt , ' dt_stiff' )

# Run solver
for i in range(n_tsteps):
    info.t0          = t[i]
    info.tf          = t[i+1]
    chparams.T_const = ch.DoubleVector(T[i] * np.ones(nx**2))
    print( 't0 = ', t[i]/lin_dt, ' dt_lin , tf = ', t[i+1]/lin_dt, ' dt_lin')
    ch.run_ch_solver(chparams,info);
