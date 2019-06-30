import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

t0 = 0
tf = 60
#statefile   = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/build/C_'
statefile   = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/swig/C_'

nx = int(np.sqrt(np.genfromtxt(statefile + '0.out' ).shape[0]))

x     = np.arange(nx)
xx,yy = np.meshgrid(x,x)

tstep = 1
fig   = plt.figure(10,figsize=(8,8))
ax    = fig.gca()

# animation function.  This is called sequentially
def animate(i):
    print(i*tstep)
    w = np.genfromtxt(statefile + str(int(i*tstep)) + '.out' )
    w = w.reshape([nx,nx]);
    ax.cla()
    contour = ax.contourf(xx,yy,w,30,vmin=-0.7,vmax=0.7)
    #contour = ax.contourf(xx,yy,w,30,vmin=0.5,vmax=1)
    return contour

anim = animation.FuncAnimation(fig, animate,
                               frames=int(tf/tstep)+1, interval=20)
anim.save('ch2d.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
#plt.show()

