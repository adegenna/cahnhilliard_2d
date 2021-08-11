from utils_plotting import *


simdata = SimData( int(input('nx : ')) , int(input('ny : ')) , int(input('nfiles : ')) )

plot_opt = input('enter plot choice from [ snapshot , video ] : ')

if plot_opt == 'video':
    make_simdata_video( simdata )

elif plot_opt == 'snapshot':
    plot_state_snapshot( input('snapshot number : ') , simdata )
    plt.show()