import numpy as np
import matplotlib.pyplot as plt

statefile   = '/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/output/state.csv'
M           = 128
N           = 128

# State plot, colored by time
C     = np.genfromtxt(statefile)
samps = C.shape[0]
plt.figure(1)
plt.ion()
for i in range(samps):
    Ci = np.reshape(C[i],[M,N])
    plt.clf()
    plt.contourf(Ci,30,vmin=-1,vmax=1)
    plt.title(i)
    plt.pause(0.01)

plt.ioff()
plt.show()

