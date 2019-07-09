import numpy as np
import matplotlib.pyplot as plt
from scipy import misc
import cahnhilliard as ch

chparams          = ch.CHparamsVector();
info              = ch.SimInfo();

# *********** INPUTS ***********
info.t0       = 0.0;
info.nx       = 128;
info.ny       = 128;
info.dx       = 1./info.nx;
info.dy       = 1./info.ny;
info.bc       = 'periodic'
info.BC_dirichlet_ch = 1.0

nx_ref        = 128;
dx_ref        = 1./nx_ref

eps_2        = 0.01**2
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

n_dt = 300
# ******************************

# Define timescales
biharm_dt         = (info.dx**4) / np.max(chparams.eps_2)
diff_dt           = (info.dx**2) / np.max( [np.max(chparams.u) , np.max(chparams.b)] )
lin_dt            = 1.0 / np.max(chparams.sigma)
biharm_dt_ref     = (dx_ref**4) / np.max(chparams.eps_2)
diff_dt_ref       = (dx_ref**2) / np.max( [np.max(chparams.u) , np.max(chparams.b)] )
lin_dt_ref        = 1.0 / np.max(chparams.sigma)

# Reset from saved state
n_tsteps          = 25
info.t0           = 0
t                 = np.linspace(info.t0 , info.t0 + n_dt * biharm_dt_ref , n_tsteps+1)
chparams.dt_check = t[1]-t[0]

# Run solver
print( 'Biharmonic timescale dt_biharm = ' , biharm_dt )
print( 'Diffusion timescale dt_diff = ' , diff_dt , ' = ' , diff_dt/biharm_dt , ' dt_biharm')
print( 'Linear timescale dt_lin = ' , lin_dt , ' = ' , lin_dt/biharm_dt , ' dt_biharm')
print( 'Sampling interval = ' , chparams.dt_check / biharm_dt , ' dt_biharm' )

for i in range(n_tsteps):
    info.t0        = t[i]
    info.tf        = t[i+1]
    rhs            = ch.CahnHilliard2DRHS(chparams , info)
    print( 't0 = ', t[i]/biharm_dt, ' dt_biharm , tf = ', t[i+1]/biharm_dt, ' dt_biharm')
    ch.run_ch_vector_nonthermal(chparams,info,rhs);
