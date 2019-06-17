import numpy as np
import matplotlib.pyplot as plt
from scipy import misc
import cahnhilliard as ch

gam0              = 0.01**2

chparams          = ch.CHparamsVector();
info              = ch.SimInfo();

info.t0       = 0.0;
info.nx       = 128;
info.dx       = 1./chparams.nx;

D        = 1.0
gamma    = gam0
b        = gamma / info.dx**2
u        = gamma / info.dx**2
alpha    = gamma * D / info.dx**4 / 200.
phi_star = 0.0;
sigma    = 0.0;

nx     = int(info.nx)
xx,yy  = np.meshgrid( np.ones(nx), np.ones(nx) )
phi_xy = xx
phi_xy[nx//2-20:nx//2-20+41 , nx//2-64:nx//2+63] = misc.imread('../../data/field_41_127.png')[:,:,-1] / 255.0
phi_xy = np.flipud(phi_xy)
#phi_xy[nx//4:3*nx//4 , nx//4:3*nx//4] = -1

chparams.D        = ch.DoubleVector(D      * np.ones(nx**2))
chparams.gamma    = ch.DoubleVector(gamma  * np.ones(nx**2))
chparams.b        = ch.DoubleVector(b      * np.ones(nx**2))
chparams.u        = ch.DoubleVector(u      * np.ones(nx**2))
chparams.alpha    = ch.DoubleVector(alpha  * np.ones(nx**2))
chparams.phi_star = ch.DoubleVector(phi_xy.ravel())
chparams.sigma    = sigma


biharm_dt         = (info.dx**4) / np.max(chparams.D) / np.max(chparams.gamma)
diff_dt           = (info.dx**2) / np.max(chparams.D) / np.max( [np.max(chparams.u) , np.max(chparams.b)] )
lin_dt            = 1.0 / np.max(chparams.alpha)
n_tsteps          = 10
t                 = np.linspace(0,300*biharm_dt,n_tsteps+1)
chparams.dt_check = t[1]-t[0]

print( 'Biharmonic timescale dt_biharm = ' , biharm_dt )
print( 'Diffusion timescale dt_diff = ' , diff_dt , ' = ' , diff_dt/biharm_dt , ' dt_biharm')
print( 'Linear timescale dt_lin = ' , lin_dt , ' = ' , lin_dt/biharm_dt , ' dt_biharm')

print( 'Sampling interval = ' , chparams.dt_check / biharm_dt , ' dt_biharm' )

for i in range(n_tsteps):
    chparams.t0        = t[i]
    chparams.tf        = t[i+1]
    print( 't0 = ', t[i]/biharm_dt, ' dt_biharm , tf = ', t[i+1]/biharm_dt, ' dt_biharm')
    rhs_scalar = ch.CahnHilliard2DRHS_Scalar(chparams,info)
    ch.run_ch_solver(chparams,info,rhs_scalar);
