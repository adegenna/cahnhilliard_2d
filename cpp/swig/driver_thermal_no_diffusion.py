import numpy as np
import matplotlib.pyplot as plt
from scipy import misc
import cahnhilliard as ch

chparams          = ch.CHparamsVector();
info              = ch.SimInfo();

# *********** INPUTS ***********
nx_ref        = 128;
dx_ref        = 1./nx_ref

info.t0       = 0.0;
info.nx       = nx_ref
info.ny       = nx_ref
info.dx       = 1./info.nx;
info.dy       = 1./info.ny;
info.bc       = 'neumann'
info.BC_dirichlet_ch = 1.0

eps_2        = 0.05**2
sigma        = eps_2 / dx_ref**4 / 200.
b            = eps_2 / dx_ref**2
u            = eps_2 / dx_ref**2
m            = 0.0;
sigma_noise  = 0.0;
DT           = 0.1*eps_2 / dx_ref**2

# Set up grid for spatial-field quantities
nx                = int(info.nx)
xx,yy             = np.meshgrid( np.arange(0,1,1/info.nx), np.arange(0,1,1/info.nx) )

chparams.eps_2        = ch.DoubleVector(eps_2  * np.ones(nx**2))
chparams.b            = ch.DoubleVector(b      * np.ones(nx**2))
chparams.u            = ch.DoubleVector(u      * np.ones(nx**2))
chparams.sigma        = ch.DoubleVector(sigma  * np.ones(nx**2))
chparams.m            = ch.DoubleVector(m      * np.ones(nx**2))
chparams.DT           = ch.DoubleVector(DT     * np.ones(nx**2))
chparams.sigma_noise  = sigma_noise
chparams.eps2_min     = 1./10*eps_2
chparams.eps2_max     =     1*eps_2
chparams.sigma_min    = 1./10*sigma
chparams.sigma_max    =     1*sigma
chparams.T_min        = 0.0
chparams.T_max        = 1.0
chparams.T_const      = ch.DoubleVector(0.  * np.ones(nx**2))

n_dt = 20000
# ******************************

# Define timescales and temporal checkpoints 
biharm_dt         = (info.dx**4) / np.max(chparams.eps_2)
diff_dt           = (info.dx**2) / np.max( [np.max(chparams.u) , np.max(chparams.b)] )
lin_dt            = 1.0 / np.max(chparams.sigma)
biharm_dt_ref     = (dx_ref**4) / np.max(chparams.eps_2)
diff_dt_ref       = (dx_ref**2) / np.max( [np.max(chparams.u) , np.max(chparams.b)] )
lin_dt_ref        = 1.0 / np.max(chparams.sigma)
n_tsteps          = 200
info.t0           = 0
t                 = np.linspace(info.t0 , info.t0 + n_dt * biharm_dt_ref , n_tsteps+1)

# Define control profile for temperature
T                = np.ones(n_tsteps)
T[0:n_tsteps//2] = np.linspace(chparams.T_min , chparams.T_max , n_tsteps//2)
T[n_tsteps//2:]  = np.linspace(chparams.T_max , chparams.T_min , n_tsteps//2)

# Run solver
print( 'Biharmonic timescale dt_biharm = ' , biharm_dt )
print( 'Diffusion timescale dt_diff = ' , diff_dt , ' = ' , diff_dt/biharm_dt , ' dt_biharm')
print( 'Linear timescale dt_lin = ' , lin_dt , ' = ' , lin_dt/biharm_dt , ' dt_biharm')
print( 'Sampling interval = ' , (t[1]-t[0]) / biharm_dt , ' dt_biharm' )

for i in range(n_tsteps):
    info.t0          = t[i]
    info.tf          = t[i+1]
    chparams.T_const = ch.DoubleVector( T[i] * np.ones(nx**2) )
    rhs              = ch.CahnHilliard2DRHS_thermal_nodiffusion(chparams , info)
    print( 't0 = ', t[i]/biharm_dt, ' dt_biharm , tf = ', t[i+1]/biharm_dt, ' dt_biharm')
    ch.run_ch_vector_thermal_nodiffusion(chparams,info,rhs);
