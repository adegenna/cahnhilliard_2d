import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def read_swig_soln_file( statefile , tstep , i ):
    
    filenametotal = statefile + str(i) + '.out'
    print( filenametotal )
    return np.genfromtxt( filenametotal )

t0 = 0.00
tf = 1.00
statefile   = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/swig/C_'
nx = 64
ny = 64

x     = np.arange(nx)
y     = np.arange(ny)
xx,yy = np.meshgrid(x,y)
xx = xx.T; yy = yy.T

tstep = 0.01
fig   = plt.figure(10,figsize=(8,8))
ax    = fig.gca()

# animation function.  This is called sequentially

frames = []

for i in range(int((tf-t0)/tstep + 1)):
    w = read_swig_soln_file( statefile , tstep , i )
    w = w.reshape([nx,ny],order='C')

    frames.append( ( ax.contourf( xx,yy,w,30,vmin=-0.4,vmax=0.4, cmap='afmhot') ).collections )

ani = animation.ArtistAnimation(fig, frames, interval=20 , blit=True )
ani.save('ch2d.mp4' )
