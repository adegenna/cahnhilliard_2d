#include "rhs_implicit.h"
#include "boundary_conditions.h"


PetscErrorCode FormIFunction(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx) {

  // Computes residual F = Udot - RHSFunction

  PetscErrorCode ierr;
  AppCtx         *user=(AppCtx*)ctx;
  DM             da   = (DM)user->da;
  PetscInt       i,j,k,Mx,My,Mz,xs,ys,zs,xm,ym,zm;
  PetscScalar      hx,hy,hz,sx,sy,sz;
  PetscScalar    u,**uarray,**f,**udot, **eps_2_array, **sigma_array;
  Vec            localU, local_eps_2, local_sigma;
  PetscScalar l_i,l_ip1,l_im1,l_jp1,l_jm1;
  PetscScalar dxx,dyy,dzz,rhs_ijk;
  PetscScalar q_im1,q_ip1,q_jm1,q_jp1,q_km1,q_kp1,q_0;
  PetscScalar sigma,m,eps_2;

  PetscFunctionBeginUser;
  
  DMGetLocalVector(da,&localU);
  DMGetLocalVector(da,&local_eps_2);
  DMGetLocalVector(da,&local_sigma);
  
  DMDAGetInfo( da ,
               PETSC_IGNORE,
               &Mx , &My , &Mz ,
               PETSC_IGNORE , PETSC_IGNORE , PETSC_IGNORE ,
               PETSC_IGNORE ,
               PETSC_IGNORE ,
               PETSC_IGNORE , PETSC_IGNORE , PETSC_IGNORE ,
               PETSC_IGNORE );

  // NOTE: these CH eqns are dimensionless with domain length scale = 1. Physical domain size shows up in L_omega.
  hx = 1.0/(PetscScalar)(Mx-1); sx = 1.0/(hx*hx);
  hy = 1.0/(PetscScalar)(My-1); sy = 1.0/(hy*hy);
  hz = 1.0/(PetscScalar)(Mz-1); sz = 1.0/(hz*hz);
  
  DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU);
  DMGlobalToLocalEnd(da,U,INSERT_VALUES,localU);
  DMGlobalToLocalBegin(da,user->eps_2,INSERT_VALUES,local_eps_2);
  DMGlobalToLocalEnd(da,user->eps_2,INSERT_VALUES,local_eps_2);
  DMGlobalToLocalBegin(da,user->sigma,INSERT_VALUES,local_sigma);
  DMGlobalToLocalEnd(da,user->sigma,INSERT_VALUES,local_sigma);
  
  DMDAVecGetArrayRead(da,localU,&uarray);
  DMDAVecGetArrayRead(da,local_eps_2,&eps_2_array);
  DMDAVecGetArrayRead(da,local_sigma,&sigma_array);
  DMDAVecGetArray(da,F,&f);
  DMDAVecGetArray(da,Udot,&udot);

  DMDAGetCorners( da ,
                  &xs , &ys , &zs ,
                  &xm , &ym , &zm );
  
  /* Compute function over the locally owned part of the grid */
  for ( k = zs ; k < zs+ym ; k++ ) {
    for ( j = ys ; j < ys+ym ; j++ ) {
      for ( i = xs ; i < xs+xm ; i++ ) {

        set_boundary_ghost_nodes( user , uarray , Mx , My , Mz , i , j , k );
	
        // dc/dt = laplacian( c^3 - c ) - laplacian( eps_2 * laplacian(c) ) - sigma*(c - m)
        sigma = sigma_array[k][j][i];
        m     = user->m;
      
        // Term: laplacian( c^3 - c )
        l_i     = uarray[k][j][i]   * uarray[k][j][i]   * uarray[k][j][i]   - uarray[k][j][i];
        l_im1   = uarray[k][j][i-1] * uarray[k][j][i-1] * uarray[k][j][i-1] - uarray[k][j][i-1];
        l_ip1   = uarray[k][j][i+1] * uarray[k][j][i+1] * uarray[k][j][i+1] - uarray[k][j][i+1];
        l_jm1   = uarray[k][j-1][i] * uarray[k][j-1][i] * uarray[k][j-1][i] - uarray[k][j-1][i];
        l_jp1   = uarray[k][j+1][i] * uarray[k][j+1][i] * uarray[k][j+1][i] - uarray[k][j+1][i];
        l_km1   = uarray[k-1][j][i] * uarray[k-1][j][i] * uarray[k-1][j][i] - uarray[k-1][j][i];
        l_kp1   = uarray[k+1][j][i] * uarray[k+1][j][i] * uarray[k+1][j][i] - uarray[k+1][j][i];
        
        dxx     = sx * ( l_ip1 + l_im1 - 2.0 * l_i );
        dyy     = sy * ( l_jp1 + l_jm1 - 2.0 * l_i );
        dzz     = sz * ( l_kp1 + l_km1 - 2.0 * l_i );
      
        rhs_ijk  = dxx + dyy + dzz;

        // Term: laplacian( -eps_2 * laplacian( c ) ) = laplacian( q )
        q_im1   = ( -eps2[k][j][i-1] ) * ( sx * ( uarray[k][j][i-2]   + uarray[k][j][i]     - 2.0*uarray[k][j][i-1] ) + sy * ( uarray[k][j-1][i-1] + uarray[k][j+1][i-1] - 2.0*uarray[k][j][i-1] ) + sz * ( uarray[k-1][j][i-1] + uarray[k+1][j][i-1] - 2.0*uarray[k][j][i-1]  ) );
        q_ip1   = ( -eps2[k][j][i+1] ) * ( sx * ( uarray[k][j][i+2]   + uarray[k][j][i]     - 2.0*uarray[k][j][i+1] ) + sy * ( uarray[k][j-1][i+1] + uarray[k][j+1][i+1] - 2.0*uarray[k][j][i+1] ) + sz * ( uarray[k-1][j][i+1] + uarray[k+1][j][i+1] - 2.0*uarray[k][j][i+1]  ) );
        q_jm1   = ( -eps2[k][j-1][i] ) * ( sx * ( uarray[k][j-1][i-1] + uarray[k][j-1][i+1] - 2.0*uarray[k][j-1][i] ) + sy * ( uarray[k][j-2][i]   + uarray[k][j][i]     - 2.0*uarray[k][j-1][i] ) + sz * ( uarray[k+1][j-1][i] + uarray[k-1][j-1][i] - 2.0*uarray[k][j-1][i]  ) );
        q_jp1   = ( -eps2[k][j+1][i] ) * ( sx * ( uarray[k][j+1][i-1] + uarray[k][j+1][i+1] - 2.0*uarray[k][j+1][i] ) + sy * ( uarray[k][j][i]     + uarray[k][j+2][i]   - 2.0*uarray[k][j+1][i] ) + sz * ( uarray[k-1][j+1][i] + uarray[k+1][j+1][i] - 2.0*uarray[k][j+1][i]  ) );
        q_km1   = ( -eps2[k-1][j][i] ) * ( sx * ( uarray[k-1][j][i-1] + uarray[k-1][j][i+1] - 2.0*uarray[k-1][j][i] ) + sy * ( uarray[k-1][j-1][i] + uarray[k-1][j+1][i] - 2.0*uarray[k-1][j][i] ) + sz * ( uarray[k-2][j][i]   + uarray[k][j][i]     - 2.0*uarray[k-1][j][i]  ) );
        q_kp1   = ( -eps2[k+1][j][i] ) * ( sx * ( uarray[k+1][j][i-1] + uarray[k+1][j][i+1] - 2.0*uarray[k+1][j][i] ) + sy * ( uarray[k+1][j-1][i] + uarray[k+1][j+1][i] - 2.0*uarray[k+1][j][i] ) + sz * ( uarray[k][j][i]     + uarray[k+2][j][i]   - 2.0*uarray[k+1][j][i]  ) );
        q_0     = ( -eps2[k][j][i]   ) * ( sx * ( uarray[k][j][i-1]   + uarray[k][j][i+1]   - 2.0*uarray[k][j][i] )   + sy * ( uarray[k][j-1][i]   + uarray[k][j+1][i]   - 2.0*uarray[k][j][i] )   + sz * ( uarray[k-1][j][i]   + uarray[k+1][j][i]   - 2.0*uarray[k][j][i]    ) );

        dxx     = sx * ( q_im1 + q_ip1 - 2.0*q_0 );
        dyy     = sy * ( q_jm1 + q_jp1 - 2.0*q_0 );
        dzz     = sz * ( q_km1 + q_kp1 - 2.0*q_0 );
        
        rhs_ijk += dxx + dyy + dzz; // laplacian( q )
      
        // Term: -sigma*(c - m)
        rhs_ijk += -sigma * ( uarray[k][j][i] - m );

        // Form f
        if ( user->boundary == 1 ) // Neumann: reset residuals explicitly 
          f[k][j][i] = reset_boundary_residual_values_for_neumann_bc( uarray , rhs_ijk , udot[k][j][i] , Mx , My , Mz , i , j , k );

        // else if ( user->boundary == 3 ) // Bottom dirichlet, rest Neumann
        //   f[k][j][i] = reset_boundary_residual_values_for_dirichlet_bottom_neumann_remainder_bc( uarray , rhs_ijk , udot[k][j][i] , Mx , My , Mz , i , j , k );

        // else if ( user->boundary == 4 ) // Bottom/top dirichlet, rest Neumann
        //   f[k][j][i] = reset_boundary_residual_values_for_dirichlet_topandbottom_neumann_remainder_bc( uarray , rhs_ijk , udot[k][j][i] , Mx , My , Mz , i , j , k );
      
        else // Dirichlet or periodic: just compute with ghost nodes
          f[k][j][i] = udot[k][j][i] - rhs_ijk;
      }

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
