import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct
from collections import namedtuple
import image_structure # git@github.com:adegenna/image_structure.git
from scipy.spatial import Voronoi, voronoi_plot_2d

def compute_image_structure_yager( data_file , outdir , N , M , interpolation_abscissa , str_figure='/C_100_' ):
    
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
                                                              str_figure=str_figure )
    return structure_metrics

def compute_image_structure_voronoi( data_file , outdir , N , M , filter_tolerance=-0.8 , str_figure='/C_100_' ):
    
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
                                                              str_figure=str_figure )
    return structure_metrics


if __name__ == '__main__':

    filename     = '/home/adegennaro/Projects/appmath/cahnhilliard_2d/data/mc_params/mc_942/C_20.out'
    outdir       = '/home/adegennaro/Projects/appmath/cahnhilliard_2d/visualization/'
    N            = 150
    M            = 150
    interpolation_abscissa  = np.linspace(0,2,200)
    
    # (x_fft , y_fft) are the (x,y) coords of the fft
    #qs , data1D , lm_result = compute_image_structure_yager( filename , outdir , N , M , interpolation_abscissa )

    C = ( np.genfromtxt( filename ).reshape([N,M]).T )

    qs , data1D , lm_result = \
          compute_image_structure_yager( filename , outdir , N , M , interpolation_abscissa , str_figure='/C_20_' )
    
    fig = plt.figure()
    ax  = plt.gca()

    ax.contourf( C )

    plt.figure()
    plt.plot( qs , data1D )
    
    plt.show()
