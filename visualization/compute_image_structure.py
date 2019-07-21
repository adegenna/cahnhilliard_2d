import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct
from collections import namedtuple
import image_structure # git@github.com:adegenna/image_structure.git

def compute_image_structure_yager( data_file , N , M ):
    
    # Named-tuple for handling input options
    Inputs      = namedtuple('Inputs' , 'data_file data_file_type dimensions structure_function output_file nx ny')
    
    w     = np.genfromtxt( data_file )
    Ci    = w.reshape([N,M],order='C');
    inputs             = Inputs(data_file, 'csv', 2, 'fourier_yager', './structure_metrics_2d.dat' , N , M)
    structure_analysis = image_structure.src.ImageStructure.ImageStructure( inputs )
    structure_metrics  = structure_analysis.compute_structure(plot_metrics=True, outdir='./', str_figure='/C_100_' )
    return structure_metrics

if __name__ == '__main__':

    filename     = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/swig/C_110.out'
    N            = 128
    M            = 128
    structure    = compute_image_structure_yager( filename , N , M )
    
