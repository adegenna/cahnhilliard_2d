import numpy as np
from dataclasses import dataclass
import glob

class SimData:

    def __init__( self , nx , ny , nfiles , statefile='../cpp/swig/C_' ):

        self.xx , self.yy = self.make_grid( nx , ny )
        self.nfiles = nfiles
        self.solnbase = statefile
    
    @property
    def nx( self ):
        return self.xx.shape[0]

    @property
    def ny( self ):
        return self.xx.shape[1]

    @property
    def base_dir( self ):
        return '/'.join( self.solnbase.split("/")[0:-1] )

    @property
    def work_dir( self ):
        return self.solnbase.split("/")[-2]

    @property
    def base_filename( self ):
        return self.solnbase.split("/")[-1]

    @property
    def filename_extension( self ):
        return ".out"

    def make_grid( self , nx , ny ):
        
        x     = np.arange(nx)
        y     = np.arange(ny)
        xx,yy = np.meshgrid(x,y)
        xx = xx.T; yy = yy.T

        return xx , yy

    def read_swig_soln_file( self , i  ):
        
        filenametotal = self.solnbase + str(i) + self.filename_extension

        try:
            return np.genfromtxt( filenametotal ).reshape( [self.nx,self.ny] , order = 'C' )
        except:
            return None

    def read_param( self , strFile : str , idx : int ):
        
        pi = np.genfromtxt( glob.glob( self.base_dir + '/' + strFile + '*' )[0] , delimiter=',' )[:,-1] # hacky way to read list of parameter values 

        return pi[idx]

@dataclass
class PlottingOptions:

    """
    class to use with plotter script
    stores options related to having many solution workdirs to draw from for plotting
    """

    nplot : int = 25
    ntotal : int = 1024
    plot_type : str = "mc_results" # [ snapshot , video , mc_results ]
    snapshot_num   : int = 20
    mc_workdir_name : str = "mc_"


def make_random_pytorch_dataset( idx_snapshot , mc_workdir_name , n_mc , simdata_base , batchsize ):

    """
    Returns a 4D numpy array ( batchsize , features=1 , nX , nY )
    """
    
    f_mc_dir = lambda idx : simdata_base.base_dir + "/" + mc_workdir_name + str(idx) + "/" + simdata_base.base_filename

    idx_rand = iter( np.random.choice( np.arange(n_mc)+1 , n_mc ) )

    U = []

    while len(U) < batchsize:
        
        try:
            si = SimData( simdata_base.nx , simdata_base.ny , simdata_base.nfiles , f_mc_dir(next(idx_rand)) )
            U.append( si.read_swig_soln_file( idx_snapshot ) )

        except:
            pass


    return np.expand_dims( np.array(U) , axis=1 )