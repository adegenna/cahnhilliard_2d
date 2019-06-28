import numpy as np
import matplotlib.pyplot as plt
from scipy import misc
import cahnhilliard as ch

chparams          = ch.CHparamsVector();
info              = ch.SimInfo();

# *********** INPUTS ***********
info.t0       = 0.0;
info.nx       = 128;
info.dx       = 1./info.nx;

eps_2        = 0.01**2
b            = eps_2 / info.dx**2
u            = eps_2 / info.dx**2
sigma        = eps_2 / info.dx**4 / 200.
m            = 0.0;
sigma_noise  = 1.0;

# Set up grid for spatial-field quantities
nx                = int(info.nx)
xx,yy  = np.meshgrid( np.zeros(nx), np.zeros(nx) )
phi_xy_AMD = xx.copy()
phi_xy_KRC = xx.copy()
phi_xy_AMD[nx//2-20:nx//2-20+41 , nx//2-64:nx//2+63] = misc.imread('../../data/field_41_127.png')[:,:,-1] / 255.0
phi_xy_AMD = np.flipud(phi_xy_AMD)
phi_xy_KRC[nx//2-26:nx//2-26+52 , nx//2-64:nx//2+64] = misc.imread('../../data/krc_52_128.png')[:,:,-1] / 255.0
phi_xy_KRC = np.flipud(phi_xy_KRC)

chparams.eps_2        = ch.DoubleVector(eps_2  * np.ones(nx**2))
chparams.b            = ch.DoubleVector(b      * np.ones(nx**2))
chparams.u            = ch.DoubleVector(u      * np.ones(nx**2))
chparams.sigma        = ch.DoubleVector(sigma  * np.ones(nx**2))
chparams.m            = ch.DoubleVector(phi_xy_AMD.ravel())
chparams.sigma_noise  = sigma_noise

n_dt = 600
# ******************************


# Define timescales
biharm_dt         = (info.dx**4) / np.max(chparams.eps_2)
diff_dt           = (info.dx**2) / np.max( [np.max(chparams.u) , np.max(chparams.b)] )
lin_dt            = 1.0 / np.max(chparams.sigma)
n_tsteps          = 60
t                 = np.linspace(0 , n_dt * biharm_dt , n_tsteps+1)
chparams.dt_check = t[1]-t[0]

# Run solver
print( 'Biharmonic timescale dt_biharm = ' , biharm_dt )
print( 'Diffusion timescale dt_diff = ' , diff_dt , ' = ' , diff_dt/biharm_dt , ' dt_biharm')
print( 'Linear timescale dt_lin = ' , lin_dt , ' = ' , lin_dt/biharm_dt , ' dt_biharm')
print( 'Sampling interval = ' , chparams.dt_check / biharm_dt , ' dt_biharm' )

for i in range(n_tsteps):
    info.t0        = t[i]
    info.tf        = t[i+1]
    if (i > 0.2*n_tsteps):
        chparams.m        = ch.DoubleVector(phi_xy_KRC.ravel())
        chparams.sigma    = ch.DoubleVector(10 * sigma  * np.ones(nx**2))
    print( 't0 = ', t[i]/biharm_dt, ' dt_biharm , tf = ', t[i+1]/biharm_dt, ' dt_biharm')
    ch.run_ch_vector(chparams,info);
