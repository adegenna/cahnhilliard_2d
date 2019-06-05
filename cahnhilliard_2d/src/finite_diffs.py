import numpy as np

def get_stencil(u):
    assert( len(u.shape) == 2)
    nx,ny = u.shape

    ul    = -1*np.ones_like(u)
    ur    = -1*np.ones_like(u)
    ut    = -1*np.ones_like(u)
    ub    = -1*np.ones_like(u)
    uc    = -1*np.ones_like(u)

    I       = np.arange(1,nx-1).reshape(nx-2,1)
    J       = np.arange(1,ny-1).reshape(1,ny-2)
    
    uc[I,J] = u[I   , J]
    ul[I,J] = u[I   , J-1]
    ur[I,J] = u[I   , J+1]
    ut[I,J] = u[I+1 , J]
    ub[I,J] = u[I-1 , J]
    
    return dict(zip( ['l','r','t','b','c'] , [ul,ur,ut,ub,uc] ))
    
def lapl(u,dx=1.):
    assert( len(u.shape) == 2)
    nx,ny    = u.shape
    u_lapl   = np.zeros_like(u)
    u_lapl[1:nx-1 , 1:ny-1]   = \
               ( -4 * u[1:nx-1   , 1:ny-1] + \
                      u[0:nx-2   , 1:ny-1] + \
                      u[2:nx     , 1:ny-1] + \
                      u[1:nx-1   , 0:ny-2] + \
                      u[1:nx-1   , 2:ny  ] ) / dx**2
    return u_lapl

def biharm(u,dx=1.):
    assert( len(u.shape) == 2)
    nx,ny     = u.shape
    u_biharm  = np.zeros_like(u)
    u_biharm[2:nx-2 , 2:ny-2]   = \
               ( 12 * u[2:nx-2   , 2:ny-2] + \
                 1  * u[0:nx-4   , 2:ny-2] + \
                 -4 * u[1:nx-3   , 2:ny-2] + \
                 -4 * u[3:nx-1   , 2:ny-2] + \
                  1 * u[4:nx     , 2:ny-2] + \
                  1 * u[2:nx-2   , 0:ny-4] + \
                 -4 * u[2:nx-2   , 1:ny-3] + \
                 -4 * u[2:nx-2   , 3:ny-1] + \
                  1 * u[2:nx-2   , 4:ny  ] ) / dx**4
    return u_biharm

def apply_neumann_bc_d1_and_d3(u,dx=1.):
    assert( len(u.shape) == 2)
    nx,ny       = u.shape
    u[1      , 2:ny-2] = u[3      , 2:ny-2]
    u[-2     , 2:ny-2] = u[-4     , 2:ny-2]
    u[2:nx-2 , 1     ] = u[2:nx-2 , 3]
    u[2:nx-2 , -2    ] = u[2:nx-2 , -4]

    u[0      , 2:ny-2] = 2*u[1,2:ny-2]  - 2*u[3,2:ny-2]  + u[4,2:ny-2]
    u[-1     , 2:ny-2] =   u[-5,2:ny-2] - 2*u[-4,2:ny-2] + 2*u[-2,2:ny-2]
    u[2:nx-2 , 0     ] = 2*u[2:nx-2,1]  - 2*u[2:nx-2,3]  + u[2:nx-2,4]
    u[2:nx-2 , -1    ] =   u[2:nx-2,-5] - 2*u[2:nx-2,-4] + 2*u[2:nx-2,-2]
    
    return u

def apply_periodic_bc(u):
    assert( len(u.shape) == 2)
    nx,ny       = u.shape
    u[0]        = u[-4]
    u[1]        = u[-3]
    u[-2]       = u[2]
    u[-1]       = u[3]
        
    u[:,0]      = u[:,-4]
    u[:,1]      = u[:,-3]
    u[:,-2]     = u[:,2]
    u[:,-1]     = u[:,3]

    return u




if __name__ == "__main__":
    
    import matplotlib.pyplot as plt

    nx = 10; ny = 9;
    x  = np.arange(nx)
    y  = np.arange(ny)+10
    xx,yy = np.meshgrid(x,y)
    zz    = xx

    print(xx)
    print(yy)

    sx = get_stencil(xx)
    sy = get_stencil(yy)

    print(sx['c'].shape)
    for i in range(ny):
        for j in range(nx):
            print('center = (', sx['c'][i,j] , ',' , sy['c'][i,j] , ')')
            print('left = ('  , sx['l'][i,j] , ',' , sy['l'][i,j] , ')')
            print('right = (' , sx['r'][i,j] , ',' , sy['r'][i,j] , ')')
            print('top = ('   , sx['t'][i,j] , ',' , sy['t'][i,j] , ')')
            print('bottom = (', sx['b'][i,j] , ',' , sy['b'][i,j] , ')')
            print('\n')

    [ print(si , '\n', sx[si]) for si in sx ]
    [ print(si , '\n', sy[si]) for si in sy ]

    zz = np.zeros_like(xx)
    zz[2:-2,2:-2] = xx[2:-2,2:-2]**4 + yy[2:-2,2:-2]**2
    zz = apply_neumann_bc_d1_and_d3(zz)
    
    print( '\n u \n' , zz , '\n lapl(u) \n' , lapl(zz) )
    print( '\n analytic lapl(u) \n' , 12*xx**2 + 2 )

    print( '\n u \n' , zz , '\n biharm(u) \n' , biharm(zz) )
    print( '\n analytic biharm(u) \n' , 24 )
