import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct
from collections import namedtuple
import image_structure # git@github.com:adegenna/image_structure.git
from scipy.spatial import Voronoi, voronoi_plot_2d

def compute_image_structure_yager( data_file , outdir , N , M , interpolation_abscissa ):
    
    # Named-tuple for handling input options
    Inputs      = namedtuple('Inputs' , 'data_file data_file_type dimensions structure_function output_file nx ny interpolation_abscissa')
    
    w     = np.genfromtxt( data_file )
    Ci    = w.reshape([N,M],order='C');
    inputs             = Inputs(data_file, 'csv', 2, 'fourier_yager_full',
                                outdir + 'structure_metrics_2d.dat' ,
                                N , M , interpolation_abscissa )
    structure_analysis = image_structure.src.ImageStructure.ImageStructure( inputs )
    structure_metrics  = structure_analysis.compute_structure(plot_metrics=True,
                                                              outdir=outdir,
                                                              str_figure='/C_100_' )
    return structure_metrics

def compute_image_structure_voronoi( data_file , outdir , N , M , filter_tolerance=-0.8 ):
    
    # Named-tuple for handling input options
    Inputs      = namedtuple('Inputs' , 'data_file data_file_type dimensions structure_function output_file nx ny filter_tolerance')
    
    w     = np.genfromtxt( data_file )
    Ci    = w.reshape([N,M],order='C');
    inputs             = Inputs(data_file, 'csv', 2, 'voronoi',
                                outdir + 'structure_metrics_2d.dat' ,
                                N , M , filter_tolerance )
    structure_analysis = image_structure.src.ImageStructure.ImageStructure( inputs )
    structure_metrics  = structure_analysis.compute_structure(plot_metrics=True,
                                                              outdir=outdir,
                                                              str_figure='/C_100_' )
    return structure_metrics


if __name__ == '__main__':

    filename     = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/data/hcp/C_110.out'
    outdir       = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/data/hcp/'
    #filename     = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/data/mc_sinusoids/mc_11/C_100.out'
    #outdir       = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/data/mc_sinusoids/mc_11/'
    N            = 128
    M            = 128
    interpolation_abscissa  = np.linspace(0,2,200)
    
    # (x_fft , y_fft) are the (x,y) coords of the fft
    #qs , data1D , lm_result = compute_image_structure_yager( filename , outdir , N , M , interpolation_abscissa )

    C = ( np.genfromtxt( filename ).reshape([128,128]).T )
    vertices_internal , centers_internal , d_cv , vor = \
          compute_image_structure_voronoi( filename , outdir , N , M , filter_tolerance=-0.8 )
    fig = plt.figure()
    ax  = plt.gca()
    ax.contourf( C )
    voronoi_plot_2d(vor , ax)

    nvert = [ len(ri) for ri in vor.regions ]
    
    for vi in vertices_internal:
        plt.plot( vi[:,0]  , vi[:,1] , 'ro' )

    n_6     = len([vi for vi in vertices_internal if len(vi) == 6])
    n_not_6 = len([vi for vi in vertices_internal if len(vi) != 6])
    idx6         = [idxi for idxi in range(len(vertices_internal)) if len(vertices_internal[idxi])==6]
    idx_not_6    = [idxi for idxi in range(len(vertices_internal)) if len(vertices_internal[idxi])!=6]

    plt.figure()
    plt.hist( d_cv[idx6] )
    plt.hist( d_cv[idx_not_6] )
    
    print( n_not_6 )
    print( n_6 )
    print( n_not_6 / ( n_not_6 + n_6 ) )
    
    plt.show()
