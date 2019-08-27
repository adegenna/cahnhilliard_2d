import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

t0 = 0.005
tf = 0.09
#statefile   = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/build/C_'
#statefile   = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/swig/C_'
#statefile   = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/data/mcruns/mc_43/C_'
statefile   = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/implicit_petsc/examples/c_'

nx = 64
ny = 64

x     = np.arange(nx)
y     = np.arange(ny)
xx,yy = np.meshgrid(x,y)
xx = xx.T; yy = yy.T

tstep = 0.005
fig   = plt.figure(10,figsize=(8,8))
ax    = fig.gca()

# animation function.  This is called sequentially
def animate(i):
    print(i*tstep + t0)
    #w = np.genfromtxt(statefile + str(int(i*tstep)) + '.out' )
    timestamp = '{:0.4f}.dat'.format( (i*tstep + t0) )
    w = np.genfromtxt(statefile + timestamp , skip_header=5 )
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

