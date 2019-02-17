import numpy as np
from scipy.fftpack import dct
import matplotlib.pyplot as plt

# Procedural script for solver
def dct2d(x,inverse=False):
    t    = 2 if not inverse else 3
    temp = dct(x,type=t,norm='ortho').transpose()
    return dct(temp,type=t,norm='ortho').transpose()

# Spatial dimensions
N     = 128
M     = 128
delx  = 1.0/(M-1)
delx2 = delx**2
x     = np.linspace(0,1,M)

# Graphics
visual_update = 10;
type_update   = 10;

# Time parameters
t     = 0;
delt  = 0.00005;
ntmax = 250;

# CH parameters
epsilon = 0.01;
eps2    = epsilon**2;

# Eyre's scheme time-step parameter
a = 2;

# time marching update parameters
lam1 = delt/delx2;
lam2 = eps2*lam1/delx2;

# unscaled eigenvalues of the laplacian (nuemann bc)
L1   = np.tile( 2*np.cos(np.pi*np.arange(N)/(N-1)) - 2 , [M,1] ).T
L2   = np.tile( 2*np.cos(np.pi*np.arange(M)/(M-1)) - 2 , [N,1] )
Leig = L1 + L2

# scaled eigenvalues of stabilized CH update matrix
CHeig = np.ones((N,M)) - (a*lam1*Leig) + (lam2*Leig*Leig)

# scaled eigenvalues of the laplacian
Seig = lam1*Leig;

# random initial conditions
U     = (np.random.uniform(0,1,[N,M])-0.5)*0.01;
hat_U = dct2d(U);

# Main loop
it = 0; j=0;
t = 0.0;

plt.figure(1);
plt.ion()
while it < ntmax:
    plt.clf()
    plt.contourf(U,30,vmin=-1,vmax=1);
    plt.title(it)
    plt.pause(0.1)
    # Update the solution
    it = it+1
    t  = t+delt
    # compute the shifted nonlinear term
    fU = (U*U*U) - ((1+a)*U)
    # compute the right hand side in tranform space
    hat_rhs = hat_U + (Seig*dct2d(fU))
    # compute the updated solution in tranform space
    hat_U   = hat_rhs/CHeig
    # invert the cosine transform
    U       = dct2d(hat_U,inverse=True)

plt.ioff();
plt.show(U)






