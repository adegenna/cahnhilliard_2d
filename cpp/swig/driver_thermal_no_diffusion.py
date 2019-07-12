import numpy as np
import matplotlib.pyplot as plt
from scipy import misc
import cahnhilliard as ch

def convert_temperature_to_flory_huggins( T , Tmin , Tmax , Xmin , Xmax ):
    assert( (T > 1e-5) & (Tmin > 1e-5) )
    return (Xmax - Xmin) / (1./Tmin - 1./Tmax) * (1./T - 1./Tmax) + Xmin

def compute_dimensionless_ch_params_from_polymer_params( L_kuhn , m , X , N , L_omega ):
    m_scaled = (1 - m) / 2.
    eps_2    = L_kuhn**2 / ( 3 * m_scaled * (1 - m_scaled) * X * L_omega**2 )
    sigma    = 36 * L_omega**2 / ( m_scaled**2 * (1 - m_scaled)**2 * L_kuhn**2 * X * N**2 )
    return eps_2 , sigma

# ********* POLYMER PARAMETERS *********
Xmin = 0.055; Xmax = 0.5;
N        =      np.mean([ 200  , 2000])
L_repeat =      (10**-9) * np.mean([ 20   , 80  ]) # meters
n_repeat = 15
L_omega  = n_repeat * L_repeat
L_kuhn   =      (10**-9) * np.mean([ 0.5 , 3.0  ]) # meters
m        = 0.0
# **************************************

Tmin = 0.1; Tmax = 1;
T    = Tmin
X                      = convert_temperature_to_flory_huggins( T , Tmin , Tmax , Xmin , Xmax )
eps2_base , sigma_base = compute_dimensionless_ch_params_from_polymer_params( L_kuhn , m , X , N , L_omega )
print( "X = " + str(X) + "\nXN = " + str(X*N) + "\nepsilon = " + str(eps2_base**0.5) + "\nsigma = " + str(sigma_base) + "\nsigma/eps^2 = " + str(sigma_base/eps2_base) )


chparams          = ch.CHparamsVector();
info              = ch.SimInfo();

# *********** INPUTS ***********
nx_ref        = 128;
dx_ref        = 1./nx_ref

info.t0       = 0.0
info.nx       = 150
info.ny       = 150
info.dx       = 1./info.nx
info.dy       = 1./info.ny
info.bc       = 'neumann'

eps_2        = eps2_base
sigma        = sigma_base
b            = 1.0
u            = 1.0
sigma_noise  = 0.0

# Set up grid for spatial-field quantities
nx                = int(info.nx)
xx,yy             = np.meshgrid( np.arange(0,1,1/info.nx), np.arange(0,1,1/info.nx) )

chparams.eps_2        = ch.DoubleVector(eps_2  * np.ones(nx**2))
chparams.b            = ch.DoubleVector(b      * np.ones(nx**2))
chparams.u            = ch.DoubleVector(u      * np.ones(nx**2))
chparams.sigma        = ch.DoubleVector(sigma  * np.ones(nx**2))
chparams.m            = ch.DoubleVector(m      * np.ones(nx**2))
chparams.sigma_noise  = sigma_noise

n_dt = 1000
# ******************************

# Define timescales and temporal checkpoints 
biharm_dt         = (info.dx**4) / np.max(chparams.eps_2)
diff_dt           = (info.dx**2) / np.max( [np.max(chparams.u) , np.max(chparams.b)] )
lin_dt            = 1.0 / np.max(chparams.sigma)

n_tsteps          = 20
info.t0           = 0
stiff_dt          = np.min([ biharm_dt , diff_dt , lin_dt ])
t                 = np.linspace(info.t0 , info.t0 + n_dt * stiff_dt , n_tsteps+1)

# Define control profile for temperature
T                = np.ones(n_tsteps)
T[0:n_tsteps//2] = np.linspace(chparams.T_min , chparams.T_max , n_tsteps//2)
T[n_tsteps//2:]  = np.linspace(chparams.T_max , chparams.T_min , n_tsteps//2)

# Run solver
print( 'Biharmonic timescale dt_biharm = ' , biharm_dt )
print( 'Diffusion timescale dt_diff = ' , diff_dt , ' = ' , diff_dt/biharm_dt , ' dt_biharm')
print( 'Linear timescale dt_lin = ' , lin_dt , ' = ' , lin_dt/biharm_dt , ' dt_biharm')
print( 'Sampling interval = ' , (t[1]-t[0]) / stiff_dt , ' dt_stiff' )

for i in range(n_tsteps):
    info.t0          = t[i]
    info.tf          = t[i+1]
    chparams.T_const = ch.DoubleVector( T[i] * np.ones(nx**2) )
    rhs              = ch.CahnHilliard2DRHS_thermal_nodiffusion(chparams , info)
    print( 't0 = ', t[i]/lin_dt, ' dt_lin , tf = ', t[i+1]/lin_dt, ' dt_lin')
    ch.run_ch_vector_thermal_nodiffusion(chparams,info,rhs);

# Save control (temperature) history
np.savetxt( 'T.out' , np.vstack( [t/lin_dt , T] ) )
