import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct
from collections import namedtuple
import image_structure # git@github.com:adegenna/image_structure.git

def dct2d(x,inverse=False):
    t    = 2 if not inverse else 3
    temp = dct(x,type=t,norm='ortho').transpose()
    return dct(temp,type=t,norm='ortho').transpose()

mcsamps     = 50
mcstart     = 50
N           = 128
M           = N
Tfinal      = 100
Tsave       = 1
basedir     = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/swig/mc_'
Tmin        = 0.1
Tmax        = 1

# Named-tuple for handling input options
Inputs             = namedtuple('Inputs' , 'data_file data_file_type dimensions structure_function output_file nx ny')

# Sort runs by temperature, in increasing order
idxT = np.zeros(mcsamps)
intT = np.zeros(mcsamps)
for i in range(mcstart , mcstart + mcsamps):
    Ti              = np.genfromtxt( basedir + str(i+1) + '/T_mc_' + str(int(i+1)) )
    intT[i-mcstart] = np.mean(Ti)

idxT = np.argsort(intT) + mcstart + 1

fig   = plt.figure(1,figsize=(10,5))
fig2  = plt.figure(2,figsize=(10,5))
fig3  = plt.figure(3,figsize=(10,5))
plt.ioff()

for i in range(mcstart , mcstart + mcsamps):
    idxI  = idxT[i-mcstart]
    print( str(idxI) + ' / ' + str(mcsamps+mcstart) )
    data_file = basedir + str(idxI) + '/C_' + str(int(Tfinal)) + '.out'
    w     = np.genfromtxt( data_file )
    Ci    = w.reshape([N,M],order='C');
    Ti    = np.genfromtxt( basedir + str(idxI) + '/T_mc_' + str(int(idxI)) )
    inputs             = Inputs(data_file, 'csv', 2, 'fourier_yager', './structure_metrics_2d.dat' , N , M)
    structure_analysis = image_structure.src.ImageStructure.ImageStructure( inputs )
    structure_metrics  = structure_analysis.compute_structure(plot_metrics=True, outdir='./', str_figure='/C_' + str(int(idxI)) )
    
    ax1 = fig.add_subplot(10,5,i+1-mcstart)
    ax1.cla()
    ax1.contourf(Ci,30,vmin=-1,vmax=1)
    ax1.set_title('C_' + str(idxI))

    ax2 = fig2.add_subplot(10,5,i+1-mcstart)
    ax2.cla()
    ax2.plot(Ti)
    ax2.set_ylim([Tmin,Tmax])
    ax2.set_title('T_' + str(idxI))
    
    # ax3 = fig3.add_subplot(10,5,i+1-mcstart)
    # ax3.cla()
    # ax3.contourf(np.log10(np.abs(dct_c)),30,vmin=-2,vmax=2)
    # ax3.set_title('log10(DCT(C))')

    plt.pause(0.01)

fig.savefig(  'C.png' )
fig2.savefig( 'T.png' )
fig3.savefig( 'DCT_C.png' )

plt.show()

