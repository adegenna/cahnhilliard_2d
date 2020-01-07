import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

def read_mpi_soln_file( statefile , timestamp , n ):
    
    c = np.zeros( n )
    count = 0
    with open( statefile + timestamp ) as f:
        for line in f:
            try:
                c[count] = line
                count += 1
            except:
                pass
    
    return c

t0 = 0.00
tf = 0.32
#statefile   = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/build/C_'
#statefile   = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/swig/C_'
statefile   = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/runs/c_'
#statefile   = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/implicit_petsc/examples/c_'

nx = 64
ny = 64

x     = np.arange(nx)
y     = np.arange(ny)
xx,yy = np.meshgrid(x,y)
xx = xx.T; yy = yy.T

tstep = 0.04
fig   = plt.figure(10,figsize=(8,8))
ax    = fig.gca()

# animation function.  This is called sequentially
def animate(i):
    print(i*tstep + t0)
    #w = np.genfromtxt(statefile + str(int(i*tstep)) + '.out' )
    timestamp = '{:0.4f}.ascii'.format( (i*tstep + t0) )
    w = read_mpi_soln_file( statefile , timestamp , nx*ny )
    #w = np.genfromtxt(statefile + timestamp , skip_header=3 )
    w = w.reshape([nx,ny],order='C');
    ax.cla()
    contour = ax.contourf(xx,yy,w,30,vmin=-0.4,vmax=0.4)
    ax.set_aspect('equal')
    #contour = ax.contourf(xx,yy,w,30,vmin=0.,vmax=1)
    return contour

anim = animation.FuncAnimation(fig, animate,
                               frames=int((tf-t0)/tstep + 1), interval=20)
anim.save('ch2d.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
#plt.show()

