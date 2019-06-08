import numpy as np
import matplotlib.pyplot as plt
import cahnhilliard as ch

gam0              = 0.01**2

chparams          = ch.CHparams();
chparams.param_type = 1 # 0 for scalar params, 1 for field params
chparams.t0       = 0.0;
chparams.nx       = 128;
chparams.dx       = 1./chparams.nx;

chparams.m        = 1.0
chparams.gam      = gam0
chparams.b        = 1.0;
chparams.u        = 1.0;
chparams.alpha    = 10.0;
chparams.phi_star = 0.0;
chparams.sigma    = 0.0;

nx = int(chparams.nx)
xx,yy = np.meshgrid( np.ones(nx), np.ones(nx) )
phi_xy = xx
phi_xy[nx//4:3*nx//4 , nx//4:3*nx//4] = -1

chparams.m_xy        = ch.DoubleVector(1.0 * np.ones(nx**2))
chparams.gam_xy      = ch.DoubleVector(gam0 * np.ones(nx**2))
chparams.b_xy        = ch.DoubleVector(1.0  * np.ones(nx**2))
chparams.u_xy        = ch.DoubleVector(1.0  * np.ones(nx**2))
chparams.alpha_xy    = ch.DoubleVector(10.0 * np.ones(nx**2))
chparams.phi_star_xy = ch.DoubleVector(phi_xy.ravel())
chparams.sigma_xy    = ch.DoubleVector(0.0  * np.ones(nx**2))


stab_dt           = 0.5*(chparams.dx**4) / chparams.m / chparams.gam
n_tsteps          = 10
t                 = np.linspace(0,200*stab_dt,n_tsteps+1)
chparams.dt_check = t[1]-t[0]

print( 'Biharmonic timescale dt_biharm = ' , stab_dt )
print( 'Sampling interval = ' , chparams.dt_check / stab_dt , ' dt_biharm' )

for i in range(n_tsteps):
    chparams.t0        = t[i]
    chparams.tf        = t[i+1]
    print( 't0 = ', t[i], ' , tf = ', t[i+1]) #, ' , phi_star = ' , chparams.phi_star )
    ch.run_ch_solver(chparams);
