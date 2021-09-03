from utils_plotting import *
from dataclasses import dataclass

@dataclass
class PlottingOptions:

    nx : int = 128
    ny : int = 128
    nfiles : int = 100
    file_base_name : str = "../data/mcruns/C_"
    plot_type : str = "mc_results" # [ snapshot , video , mc_results ]
    snapshot_num   : int = 50
    mc_workdir_name : str = "mc_"
    n_mc : int = 100




def main( popt : PlottingOptions ):

    simdata = SimData( popt.nx , popt.ny , popt.nfiles , popt.file_base_name )

    if popt.plot_type == 'video':
        make_simdata_video( simdata )

    elif popt.plot_type == 'snapshot':
        plot_state_snapshot( popt.snapshot_num , simdata )
        plt.show()

    elif popt.plot_type == 'mc_results':
        plot_mc_results( popt.snapshot_num , \
                         popt.mc_workdir_name , \
                         popt.n_mc , \
                         simdata )



if __name__ == "__main__":

    main( PlottingOptions(snapshot_num=2) )
