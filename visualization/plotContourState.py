import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct

def dct2d(x,inverse=False):
    t    = 2 if not inverse else 3
    temp = dct(x,type=t,norm='ortho').transpose()
    return dct(temp,type=t,norm='ortho').transpose()

statefile   = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/output/C_'
M           = 128
N           = 128
Tfinal      = 100000
Tsave       = 5000

# Read data
samps = Tfinal/Tsave + 1
C     = np.zeros([samps,M*N])
for i in range(samps):
    timestamp = int(Tsave*i)
    C[i]      = np.genfromtxt(statefile + str(timestamp) + '.out')

# Precompute fourier cosine spectrum
dct_c  = np.zeros([samps,M,N])
for i in range(samps):
    Ci       = np.reshape(C[i],[M,N])
    dct_c[i] = dct2d(Ci)

# State plot, colored by time
fig   = plt.figure(1,figsize=(10,5))
ax1   = fig.add_subplot(221)
ax2   = fig.add_subplot(222)
ax3   = fig.add_subplot(223)
plt.ion()
for i in range(samps):
    Ci = np.reshape(C[i],[M,N])
    if (i != 0):
        Cim1 = np.reshape(C[i-1],[M,N])
        dcdt = np.linalg.norm(Ci-Cim1)
    else:
        dcdt = 0
    ax1.cla()
    ax1.contourf(Ci,30,vmin=-1,vmax=1)
    ax1.set_title(i)
    ax2.scatter(i+1,dcdt,c='b')
    ax2.set_yscale('log')
    ax2.set_ylim([1,100])
    ax2.set_title('dC/dt')
    ax3.cla()
    ax3.contourf(np.log10(np.abs(dct_c[i])),30,vmin=-2,vmax=2)
    ax3.set_title('log10(DCT(C))')
    plt.pause(0.01)

# Static plot of 3 different times
plt.ioff()
fig   = plt.figure(2,figsize=(10,5))
ax1   = fig.add_subplot(231)
ax2   = fig.add_subplot(234)
ax3   = fig.add_subplot(232)
ax4   = fig.add_subplot(235)
ax5   = fig.add_subplot(233)
ax6   = fig.add_subplot(236)
ax    = [ax1,ax2,ax3,ax4,ax5,ax6]
T     = [0,samps/2,samps-1]
tr    = 32
for i in range(len(T)):
    ax[2*i].contourf( np.reshape(C[T[i]],[M,N]) , 30, vmin=-1, vmax=1 )
    ax[2*i+1].contourf( np.log10(np.abs(dct_c[T[i]]))[0:tr,0:tr],30,vmin=-1,vmax=1 )
    ax[2*i].set_title(T[i])

plt.show()
