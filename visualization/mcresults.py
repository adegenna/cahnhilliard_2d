import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct
from collections import namedtuple
from scipy.spatial import Voronoi , voronoi_plot_2d
import image_structure # git@github.com:adegenna/image_structure.git

def dct2d(x,inverse=False):
    t    = 2 if not inverse else 3
    temp = dct(x,type=t,norm='ortho').transpose()
    return dct(temp,type=t,norm='ortho').transpose()

def compute_voronoi_lattice_metrics( data_file , N , M ):

    w     = np.genfromtxt( data_file )
    Ci    = w.reshape([N,M],order='C');
    inputs             = Inputs(data_file, 'csv', 2, 'voronoi', './structure_metrics_2d.dat' , N , M, filter_tolerance=-0.8)
    structure_analysis = image_structure.src.ImageStructure.ImageStructure( inputs )
    vertices_internal , centers_internal , vol_internal , d_cv , vor = \
                structure_analysis.compute_structure(plot_metrics=False, outdir='./', str_figure='/C_' )
    n_6           = len([vi for vi in vertices_internal if len(vi) == 6])
    n_not_6       = len([vi for vi in vertices_internal if len(vi) != 6])
    lattice_ratio = n_not_6 / ( n_not_6 + n_6 )

    return vor , lattice_ratio , vol_internal , centers_internal


mcsamps     = 50
mcstart     = 0
N           = 128
M           = N
Tfinal      = 100
Tsave       = 1
basedir     = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/data/mcruns/mc_'
Tmin        = 0.1
Tmax        = 1

# Named-tuple for handling input options
Inputs             = namedtuple('Inputs' , 'data_file data_file_type dimensions structure_function output_file nx ny filter_tolerance')

# Sort runs by temperature, in increasing order
idxT  = np.zeros(mcsamps)
intT  = np.zeros(mcsamps)
k     = np.zeros(mcsamps)
idxk  = np.zeros(mcsamps)
n_samples = 0
for i in range(mcstart , mcstart + mcsamps):
    try:
        Ti              = np.genfromtxt( basedir + str(i+1) + '/T_mc_' + str(int(i+1)) + '.out')
        data_file       = basedir + str(i+1) + '/C_' + str(int(Tfinal)) + '.out'
        intT[i-mcstart] = np.mean(Ti)
        omega           = np.abs(np.fft.fft( Ti ) )[1:len(Ti)//2]
        freq            = np.arange( 1 , len(Ti)//2 )
        k[i-mcstart]    = np.sum( omega * freq )
        n_samples      += 1
    except:
        print('skipping ' + str(i+1) )
        pass
    
idxT = np.argsort(intT) + mcstart + 1
idxk = np.argsort(k)    + mcstart + 1

fig   = plt.figure(1,figsize=(15,7))
fig2  = plt.figure(2,figsize=(15,7))
fig3,ax3  = plt.subplots(1,3,figsize=(10,5))

subplotsx = int( np.ceil(np.sqrt( n_samples//2 ) ) )
subplotsy = n_samples // subplotsx
plt.ioff()

x       = np.arange(N)
xx , yy = np.meshgrid( x , x )
xx      = xx.T
yy      = yy.T

plot_count = 0
for i in range(mcstart , mcstart + mcsamps):
    try:
        idxI  = idxT[i-mcstart]
        print( str(idxI) + ' / ' + str(mcsamps+mcstart) )
        data_file     = basedir + str(idxI) + '/C_' + str(int(Tfinal)) + '.out'
        w             = np.genfromtxt( data_file )
        Ci            = w.reshape([N,M],order='C');
        Ti            = np.genfromtxt( basedir + str(idxI) + '/T_mc_' + str(int(idxI)) + '.out')
        vor , lattice_ratio , vol_internal , centers_internal = compute_voronoi_lattice_metrics( data_file , N , M )

        print( vor.vertices)
        print( vor.regions )
        
        ax1 = fig.add_subplot(subplotsx , subplotsy , plot_count+1)
        ax1.cla()
        ax1.contourf(xx,yy,Ci,30,vmin=-1,vmax=1)
        ax1.plot( np.array(centers_internal)[:,0] , np.array(centers_internal)[:,1] , 'ro' )
        voronoi_plot_2d( vor , ax1 , show_vertices=False)
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.set_title('C_' + str(idxI) + ', L=%.2f' %lattice_ratio + ', M=%.0f' %np.mean(vol_internal) + ', V=%.0f' %np.var(vol_internal) )
        
        ax2 = fig2.add_subplot(subplotsx , subplotsy , plot_count+1)
        ax2.cla()
        ax2.plot(Ti)
        ax2.set_ylim([Tmin,Tmax])
        ax2.set_title('T_' + str(idxI))

        ax3[0].plot( intT[idxI-1]      , lattice_ratio , 'bo' )
        ax3[0].set_title( 'L vs avg(T)' )
        ax3[1].plot( intT[idxI-1]      , np.mean(vol_internal) , 'bo' )
        ax3[1].set_title( 'M vs avg(T)' )
        ax3[2].semilogy( intT[idxI-1]  , np.var(vol_internal) , 'bo' )
        ax3[2].set_title( 'log(V) vs avg(T)' )
        
        plot_count += 1
        
        plt.pause(0.01)
    except:
        print( 'skipping ' + str(idxI) )
        pass

fig.savefig(  'C.png' )
fig2.savefig( 'T.png' )
fig3.savefig( 'lattice.png' )

plt.show()

