#include "rhs_implicit.h"
#include "boundary_conditions.h"


PetscErrorCode FormIFunction(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx) {

  // Computes residual F = Udot - RHSFunction

  PetscErrorCode ierr;
  AppCtx         *user=(AppCtx*)ctx;
  DM             da   = (DM)user->da;
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscReal      hx,hy,sx,sy;
  PetscScalar    u,uxx,uyy,**uarray,**f,**udot, **eps_2_array, **sigma_array;
  Vec            localU, local_eps_2, local_sigma;
  
  PetscReal c_i,c_im2,c_im1,c_ip1,c_ip2,c_jm2,c_jm1,c_jp1,c_jp2, c_ul,c_ur,c_bl,c_br;
  PetscReal l_i,l_ip1,l_im1,l_jp1,l_jm1;
  PetscReal dxx,dyy,dxxxx,dyyyy,dxxyy,rhs_ij;
  PetscScalar q_im1,q_ip1,q_jm1,q_jp1,q_0;
  PetscReal sigma,m,eps_2;

  PetscFunctionBeginUser;
  DMGetLocalVector(da,&localU);
  DMGetLocalVector(da,&local_eps_2);
  DMGetLocalVector(da,&local_sigma);
  
  DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);

  // NOTE: these CH eqns are dimensionless with domain length scale = 1. Physical domain size shows up in L_omega.
  hx = 1.0/(PetscReal)(Mx-1); sx = 1.0/(hx*hx);
  hy = 1.0/(PetscReal)(My-1); sy = 1.0/(hy*hy);

  /*
     Scatter ghost points to local vector,using the 2-step process
        DMGlobalToLocalBegin(),DMGlobalToLocalEnd().
     By placing code between these two statements, computations can be
     done while messages are in transition.
  */
  DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU);
  DMGlobalToLocalEnd(da,U,INSERT_VALUES,localU);
  DMGlobalToLocalBegin(da,user->eps_2,INSERT_VALUES,local_eps_2);
  DMGlobalToLocalEnd(da,user->eps_2,INSERT_VALUES,local_eps_2);
  DMGlobalToLocalBegin(da,user->sigma,INSERT_VALUES,local_sigma);
  DMGlobalToLocalEnd(da,user->sigma,INSERT_VALUES,local_sigma);
  
  /* Get pointers to vector data */
  DMDAVecGetArrayRead(da,localU,&uarray);
  DMDAVecGetArrayRead(da,local_eps_2,&eps_2_array);
  DMDAVecGetArrayRead(da,local_sigma,&sigma_array);
  DMDAVecGetArray(da,F,&f);
  DMDAVecGetArray(da,Udot,&udot);

  /* Get local grid boundaries */
  DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);

  /* Compute function over the locally owned part of the grid */
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {

      set_boundary_ghost_nodes( user , uarray , Mx , My , i , j );

      c_i   = uarray[j][i];
      
      /* Boundary conditions */
      ThirteenPointStencil stencil      = get_thirteen_point_stencil( user , uarray      , Mx , My , i , j );
      ThirteenPointStencil stencil_eps2 = get_thirteen_point_stencil( user , eps_2_array , Mx , My , i , j );
      
      // dc/dt = laplacian( c^3 - c ) - eps_2*biharm(c) - sigma*(c - m)
      sigma = sigma_array[j][i];
      m     = user->m;
      
      // Term: laplacian( c^3 - c )
      l_i     = stencil.c_i*stencil.c_i*stencil.c_i       - stencil.c_i;
      l_im1   = stencil.c_im1*stencil.c_im1*stencil.c_im1 - stencil.c_im1;
      l_ip1   = stencil.c_ip1*stencil.c_ip1*stencil.c_ip1 - stencil.c_ip1;
      l_jm1   = stencil.c_jm1*stencil.c_jm1*stencil.c_jm1 - stencil.c_jm1;
      l_jp1   = stencil.c_jp1*stencil.c_jp1*stencil.c_jp1 - stencil.c_jp1;

      dxx     = sx * ( l_ip1 + l_im1 - 2.0 * l_i );
      dyy     = sy * ( l_jp1 + l_jm1 - 2.0 * l_i );
      
      rhs_ij  = dxx + dyy;

      // Term: laplacian( -eps_2 * laplacian( c ) ) = laplacian( q )
      q_im1   = ( -stencil_eps2.c_im1 ) * ( sx * ( stencil.c_im2 + stencil.c_i   - 2.0*stencil.c_im1 ) + sy * ( stencil.c_ul  + stencil.c_bl  - 2.0*stencil.c_im1 ) );
      q_ip1   = ( -stencil_eps2.c_ip1 ) * ( sx * ( stencil.c_i   + stencil.c_ip2 - 2.0*stencil.c_ip1 ) + sy * ( stencil.c_ur  + stencil.c_br  - 2.0*stencil.c_ip1 ) );
      q_jm1   = ( -stencil_eps2.c_jm1 ) * ( sx * ( stencil.c_ul  + stencil.c_ur  - 2.0*stencil.c_jm1 ) + sy * ( stencil.c_jm2 + stencil.c_i   - 2.0*stencil.c_jm1 ) );
      q_jp1   = ( -stencil_eps2.c_jp1 ) * ( sx * ( stencil.c_bl  + stencil.c_br  - 2.0*stencil.c_jp1 ) + sy * ( stencil.c_jp2 + stencil.c_i   - 2.0*stencil.c_jp1 ) );
      q_0     = ( -stencil_eps2.c_i   ) * ( sx * ( stencil.c_im1 + stencil.c_ip1 - 2.0*stencil.c_i   ) + sy * ( stencil.c_jm1 + stencil.c_jp1 - 2.0*stencil.c_i   ) );
      
      rhs_ij += sx * ( q_im1 + q_ip1 - 2.0*q_0 ) + sy * ( q_jm1 + q_jp1 - 2.0*q_0 ); // laplacian( q )
      
      // Term: -eps_2*biharm(c)
      // dxxxx = sx * sx * (stencil.c_ip2 - 4.0*stencil.c_ip1 + 6.0*stencil.c_i - 4.0*stencil.c_im1 + stencil.c_im2);
      // dyyyy = sy * sy * (stencil.c_jp2 - 4.0*stencil.c_jp1 + 6.0*stencil.c_i - 4.0*stencil.c_jm1 + stencil.c_jm2);
      // dxxyy = sx * sy * 2 * (4*stencil.c_i - 2*(stencil.c_im1 + stencil.c_ip1 + stencil.c_jm1 + stencil.c_jp1) + stencil.c_ul + stencil.c_ur + stencil.c_bl + stencil.c_br );	// mixed term 2*u_xxyy
      // rhs_ij += -eps_2 * ( dxxxx + dyyyy + dxxyy );

      // Term: -sigma*(c - m)
      rhs_ij += -sigma * ( stencil.c_i - m );

      // Form f
      if ( user->boundary == 1 ) // Neumann: reset residuals explicitly 
        f[j][i] = reset_boundary_residual_values_for_neumann_bc( uarray , rhs_ij , udot[j][i] , Mx , My , i , j );

      else if ( user->boundary == 3 ) // Bottom dirichlet, rest Neumann
        f[j][i] = reset_boundary_residual_values_for_dirichlet_bottom_neumann_remainder_bc( uarray , rhs_ij , udot[j][i] , Mx , My , i , j );

      else if ( user->boundary == 4 ) // Bottom/top dirichlet, rest Neumann
        f[j][i] = reset_boundary_residual_values_for_dirichlet_topandbottom_neumann_remainder_bc( uarray , rhs_ij , udot[j][i] , Mx , My , i , j );
      
      else // Dirichlet or periodic: just compute with ghost nodes
        f[j][i] = udot[j][i] - rhs_ij;
    }

  }
  
  /* Restore vectors */
  DMDAVecRestoreArrayRead(da,localU,&uarray);
  DMDAVecRestoreArrayRead(da,local_eps_2,&eps_2_array);
  DMDAVecRestoreArrayRead(da,local_sigma,&sigma_array);
  DMDAVecRestoreArray(da,F,&f);
  DMDAVecRestoreArray(da,Udot,&udot);
  DMRestoreLocalVector(da,&localU);
  DMRestoreLocalVector(da,&local_eps_2);
  DMRestoreLocalVector(da,&local_sigma);
  
  PetscFunctionReturn(0);
}
