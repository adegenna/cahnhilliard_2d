#include "rhs_implicit.h"
#include "boundary_conditions.h"
#include <petscdmcomposite.h>

PetscScalar*** FormLocalResidual( DMDALocalInfo *info ,
                                 PetscScalar ***uarray ,
                                 PetscScalar ***f , 
                                 PetscScalar ***udot ,
                                 PetscScalar ***rhs ,
                                 AppCtx *user ) {
  
  // Apply correction to boundary fluxes and set residuals
  for (int k=info->zs; k<info->zs+info->zm; k++) {
    for (int j=info->ys; j<info->ys+info->ym; j++) {
      for (int i=info->xs; i<info->xs+info->xm; i++) {

        f[k][j][i] = user->residualFunction( uarray , rhs[k][j][i] , udot[k][j][i] , info->mx , info->my , info->mz , i , j , k );
      
      }
    }
  }
  
  return f;

}

// PetscErrorCode FormIFunction(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx) {

//   // Computes residual F = Udot - RHSFunction

//   PetscErrorCode ierr;
//   AppCtx         *user=(AppCtx*)ctx;
//   DM             da   = (DM)user->da;
//   PetscInt       i,j,k,Mx,My,Mz,xs,ys,zs,xm,ym,zm;
//   PetscScalar      hx,hy,hz,sx,sy,sz;
//   PetscScalar    u,***uarray,***f,***udot, ***eps_2_array, ***sigma_array;
//   Vec            localU, localUdot, localF, local_eps_2, local_sigma;
//   PetscScalar l_i,l_ip1,l_im1,l_jp1,l_jm1,l_km1,l_kp1;
//   PetscScalar dxx,dyy,dzz,rhs_ijk;
//   PetscScalar q_im1,q_ip1,q_jm1,q_jp1,q_km1,q_kp1,q_0;
//   PetscScalar sigma,m,eps_2;

//   PetscFunctionBeginUser;
  
//   DMGetLocalVector(da,&localU);
//   // DMGetLocalVector(da,&localUdot);
//   // DMGetLocalVector(da,&localF);
//   DMGetLocalVector(da,&local_eps_2);
//   DMGetLocalVector(da,&local_sigma);
  
//   DMDAGetInfo( da ,
//                PETSC_IGNORE,
//                &Mx , &My , &Mz ,
//                PETSC_IGNORE , PETSC_IGNORE , PETSC_IGNORE ,
//                PETSC_IGNORE ,
//                PETSC_IGNORE ,
//                PETSC_IGNORE , PETSC_IGNORE , PETSC_IGNORE ,
//                PETSC_IGNORE );

//   // NOTE: these CH eqns are dimensionless with domain length scale = 1. Physical domain size shows up in L_omega.
//   hx = 1.0/(PetscScalar)(Mx-1); sx = 1.0/(hx*hx);
//   hy = 1.0/(PetscScalar)(My-1); sy = 1.0/(hy*hy);
//   hz = 1.0/(PetscScalar)(Mz-1); sz = 1.0/(hz*hz);
  
//   DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU);
//   DMGlobalToLocalEnd(da,U,INSERT_VALUES,localU);
//   // DMGlobalToLocalBegin(da,Udot,INSERT_VALUES,localUdot);
//   // DMGlobalToLocalEnd(da,Udot,INSERT_VALUES,localUdot);
//   // DMGlobalToLocalBegin(da,F,INSERT_VALUES,localF);
//   // DMGlobalToLocalEnd(da,F,INSERT_VALUES,localF);
//   DMGlobalToLocalBegin(da,user->eps_2,INSERT_VALUES,local_eps_2);
//   DMGlobalToLocalEnd(da,user->eps_2,INSERT_VALUES,local_eps_2);
//   DMGlobalToLocalBegin(da,user->sigma,INSERT_VALUES,local_sigma);
//   DMGlobalToLocalEnd(da,user->sigma,INSERT_VALUES,local_sigma);
  
//   DMDAVecGetArrayRead(da,localU,&uarray);
//   DMDAVecGetArray(da,Udot,&udot);
//   DMDAVecGetArray(da,F,&f);  
//   DMDAVecGetArrayRead(da,local_eps_2,&eps_2_array);
//   DMDAVecGetArrayRead(da,local_sigma,&sigma_array);

//   DMDAGetCorners( da ,
//                   &xs , &ys , &zs ,
//                   &xm , &ym , &zm );

//   // Set boundary ghost nodes for eps_2 and U
//   for ( k = zs ; k < zs+zm ; k++ ) {
//     for ( j = ys ; j < ys+ym ; j++ ) {
//       for ( i = xs ; i < xs+xm ; i++ ) {

//         set_boundary_ghost_nodes( user , uarray , Mx , My , Mz , i , j , k );
//         set_boundary_ghost_nodes_normal_extrapolation( user , eps_2_array , Mx , My , Mz , i , j , k );

//       }
//     }
//   }
        
//   /* Compute function over the locally owned part of the grid */
//   for ( k = zs ; k < zs+zm ; k++ ) {
//     for ( j = ys ; j < ys+ym ; j++ ) {
//       for ( i = xs ; i < xs+xm ; i++ ) {
	
//         // dc/dt = laplacian( c^3 - c ) - laplacian( eps_2 * laplacian(c) ) - sigma*(c - m)
//         sigma = sigma_array[k][j][i];
//         m     = user->m;
      
//         // Term: laplacian( c^3 - c )
//         l_i     = uarray[k][j][i]   * uarray[k][j][i]   * uarray[k][j][i]   - uarray[k][j][i];
//         l_im1   = uarray[k][j][i-1] * uarray[k][j][i-1] * uarray[k][j][i-1] - uarray[k][j][i-1];
//         l_ip1   = uarray[k][j][i+1] * uarray[k][j][i+1] * uarray[k][j][i+1] - uarray[k][j][i+1];
//         l_jm1   = uarray[k][j-1][i] * uarray[k][j-1][i] * uarray[k][j-1][i] - uarray[k][j-1][i];
//         l_jp1   = uarray[k][j+1][i] * uarray[k][j+1][i] * uarray[k][j+1][i] - uarray[k][j+1][i];
//         l_km1   = uarray[k-1][j][i] * uarray[k-1][j][i] * uarray[k-1][j][i] - uarray[k-1][j][i];
//         l_kp1   = uarray[k+1][j][i] * uarray[k+1][j][i] * uarray[k+1][j][i] - uarray[k+1][j][i];
        
//         dxx     = sx * ( l_ip1 + l_im1 - 2.0 * l_i );
//         dyy     = sy * ( l_jp1 + l_jm1 - 2.0 * l_i );
//         dzz     = sz * ( l_kp1 + l_km1 - 2.0 * l_i );
      
//         rhs_ijk  = dxx + dyy + dzz;

//         // Term: laplacian( -eps_2 * laplacian( c ) ) = laplacian( q )
//         q_im1   = ( -eps_2_array[k][j][i-1] ) * ( sx * ( uarray[k][j][i-2]   + uarray[k][j][i]     - 2.0*uarray[k][j][i-1] ) + sy * ( uarray[k][j-1][i-1] + uarray[k][j+1][i-1] - 2.0*uarray[k][j][i-1] ) + sz * ( uarray[k-1][j][i-1] + uarray[k+1][j][i-1] - 2.0*uarray[k][j][i-1]  ) );
//         q_ip1   = ( -eps_2_array[k][j][i+1] ) * ( sx * ( uarray[k][j][i+2]   + uarray[k][j][i]     - 2.0*uarray[k][j][i+1] ) + sy * ( uarray[k][j-1][i+1] + uarray[k][j+1][i+1] - 2.0*uarray[k][j][i+1] ) + sz * ( uarray[k-1][j][i+1] + uarray[k+1][j][i+1] - 2.0*uarray[k][j][i+1]  ) );
//         q_jm1   = ( -eps_2_array[k][j-1][i] ) * ( sx * ( uarray[k][j-1][i-1] + uarray[k][j-1][i+1] - 2.0*uarray[k][j-1][i] ) + sy * ( uarray[k][j-2][i]   + uarray[k][j][i]     - 2.0*uarray[k][j-1][i] ) + sz * ( uarray[k+1][j-1][i] + uarray[k-1][j-1][i] - 2.0*uarray[k][j-1][i]  ) );
//         q_jp1   = ( -eps_2_array[k][j+1][i] ) * ( sx * ( uarray[k][j+1][i-1] + uarray[k][j+1][i+1] - 2.0*uarray[k][j+1][i] ) + sy * ( uarray[k][j][i]     + uarray[k][j+2][i]   - 2.0*uarray[k][j+1][i] ) + sz * ( uarray[k-1][j+1][i] + uarray[k+1][j+1][i] - 2.0*uarray[k][j+1][i]  ) );
//         q_km1   = ( -eps_2_array[k-1][j][i] ) * ( sx * ( uarray[k-1][j][i-1] + uarray[k-1][j][i+1] - 2.0*uarray[k-1][j][i] ) + sy * ( uarray[k-1][j-1][i] + uarray[k-1][j+1][i] - 2.0*uarray[k-1][j][i] ) + sz * ( uarray[k-2][j][i]   + uarray[k][j][i]     - 2.0*uarray[k-1][j][i]  ) );
//         q_kp1   = ( -eps_2_array[k+1][j][i] ) * ( sx * ( uarray[k+1][j][i-1] + uarray[k+1][j][i+1] - 2.0*uarray[k+1][j][i] ) + sy * ( uarray[k+1][j-1][i] + uarray[k+1][j+1][i] - 2.0*uarray[k+1][j][i] ) + sz * ( uarray[k][j][i]     + uarray[k+2][j][i]   - 2.0*uarray[k+1][j][i]  ) );
//         q_0     = ( -eps_2_array[k][j][i]   ) * ( sx * ( uarray[k][j][i-1]   + uarray[k][j][i+1]   - 2.0*uarray[k][j][i] )   + sy * ( uarray[k][j-1][i]   + uarray[k][j+1][i]   - 2.0*uarray[k][j][i] )   + sz * ( uarray[k-1][j][i]   + uarray[k+1][j][i]   - 2.0*uarray[k][j][i]    ) );

//         dxx     = sx * ( q_im1 + q_ip1 - 2.0*q_0 );
//         dyy     = sy * ( q_jm1 + q_jp1 - 2.0*q_0 );
//         dzz     = sz * ( q_km1 + q_kp1 - 2.0*q_0 );
        
//         rhs_ijk += dxx + dyy + dzz; // laplacian( q )
      
//         // Term: -sigma*(c - m)
//         rhs_ijk += -sigma * ( uarray[k][j][i] - m );

//         // Form f
//         if ( user->boundary == 1 ) // Neumann: reset residuals explicitly 
//           f[k][j][i] = reset_boundary_residual_values_for_neumann_bc( uarray , rhs_ijk , udot[k][j][i] , Mx , My , Mz , i , j , k );

//         // else if ( user->boundary == 3 ) // Bottom dirichlet, rest Neumann
//         //   f[k][j][i] = reset_boundary_residual_values_for_dirichlet_bottom_neumann_remainder_bc( uarray , rhs_ijk , udot[k][j][i] , Mx , My , Mz , i , j , k );

//         // else if ( user->boundary == 4 ) // Bottom/top dirichlet, rest Neumann
//         //   f[k][j][i] = reset_boundary_residual_values_for_dirichlet_topandbottom_neumann_remainder_bc( uarray , rhs_ijk , udot[k][j][i] , Mx , My , Mz , i , j , k );
      
//         else // Dirichlet or periodic: just compute with ghost nodes
//           f[k][j][i] = udot[k][j][i] - rhs_ijk;
//       }

//     }
//   }
//   /* Restore vectors */
//   DMDAVecRestoreArrayRead(da,localU,&uarray);
//   DMDAVecRestoreArray(da,Udot,&udot);
//   DMDAVecRestoreArray(da,F,&f);
//   DMDAVecRestoreArrayRead(da,local_eps_2,&eps_2_array);
//   DMDAVecRestoreArrayRead(da,local_sigma,&sigma_array);

//   // DMLocalToGlobalBegin( da , localF , INSERT_VALUES , F );
//   // DMLocalToGlobalEnd(   da , localF , INSERT_VALUES , F );
  
//   DMRestoreLocalVector(da,&localU);
//   // DMRestoreLocalVector(da,&localUdot);
//   // DMRestoreLocalVector(da,&localF);
//   DMRestoreLocalVector(da,&local_eps_2);
//   DMRestoreLocalVector(da,&local_sigma);
  
//   PetscFunctionReturn(0);

// }


PetscScalar*** FormLocalImplicitResidualTEST( DMDALocalInfo *info ,
					      PetscScalar ***uarray ,
					      PetscScalar ***f , 
					      PetscScalar ***udot ,
					      PetscScalar ***rhs ,
					      AppCtx *user ) {

  PetscScalar rhs_ijk;

  for ( int k = info->zs ; k < info->zs+info->zm ; k++ ) {
    for ( int j = info->ys ; j < info->ys+info->ym ; j++ ) {
      for ( int i = info->xs ; i < info->xs+info->xm ; i++ ) {
	rhs_ijk    = -1.0 * uarray[k][j][i];	
	f[k][j][i] = udot[k][j][i] - rhs_ijk;
      }
    }
  }
  
  return f;

}

PetscScalar*** FormLocalRHSTEST( DMDALocalInfo *info ,
				PetscScalar ***uarray ,
				PetscScalar ***rhs , 
				AppCtx *user ) {

  for ( int k = info->zs ; k < info->zs+info->zm ; k++ ) {
    for ( int j = info->ys ; j < info->ys+info->ym ; j++ ) {
      for ( int i = info->xs ; i < info->xs+info->xm ; i++ ) {
	rhs[k][j][i]    = -1.0 * uarray[k][j][i];	
      }
    }
  }
  
  return rhs;

}


// PetscErrorCode FormIFunctionTEST(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx) {

//   // Computes residual F = Udot - RHSFunction

//   PetscErrorCode ierr;
//   AppCtx         *user=(AppCtx*)ctx;
//   DM             da   = (DM)user->da;
//   PetscInt       i,j,k,Mx,My,Mz,xs,ys,zs,xm,ym,zm;
//   PetscScalar      hx,hy,hz,sx,sy,sz;
//   PetscScalar    u,***uarray,***f,***udot, ***eps_2_array, ***sigma_array;
//   Vec            localU, localUdot, localF, local_eps_2, local_sigma;
//   PetscScalar l_i,l_ip1,l_im1,l_jp1,l_jm1,l_km1,l_kp1;
//   PetscScalar dxx,dyy,dzz,rhs_ijk;
//   PetscScalar q_im1,q_ip1,q_jm1,q_jp1,q_km1,q_kp1,q_0;
//   PetscScalar sigma,m,eps_2;

//   PetscFunctionBeginUser;
  
//   DMGetLocalVector(da,&localU);
//   // DMGetLocalVector(da,&localUdot);
//   // DMGetLocalVector(da,&localF);
  
//   DMDAGetInfo( da ,
//                PETSC_IGNORE,
//                &Mx , &My , &Mz ,
//                PETSC_IGNORE , PETSC_IGNORE , PETSC_IGNORE ,
//                PETSC_IGNORE ,
//                PETSC_IGNORE ,
//                PETSC_IGNORE , PETSC_IGNORE , PETSC_IGNORE ,
//                PETSC_IGNORE );
  
//   DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU);
//   DMGlobalToLocalEnd(da,U,INSERT_VALUES,localU);
//   // DMGlobalToLocalBegin(da,Udot,INSERT_VALUES,localUdot);
//   // DMGlobalToLocalEnd(da,Udot,INSERT_VALUES,localUdot);
//   // DMGlobalToLocalBegin(da,F,INSERT_VALUES,localF);
//   // DMGlobalToLocalEnd(da,F,INSERT_VALUES,localF);
  
//   DMDAVecGetArrayRead(da,localU,&uarray);
//   DMDAVecGetArray(da,Udot,&udot);
//   DMDAVecGetArray(da,F,&f);  
  
//   DMDAGetCorners( da ,
//                   &xs , &ys , &zs ,
//                   &xm , &ym , &zm );
  
//   /* Compute function over the locally owned part of the grid */
//   for ( k = zs ; k < zs+zm ; k++ ) {
//     for ( j = ys ; j < ys+ym ; j++ ) {
//       for ( i = xs ; i < xs+xm ; i++ ) {
// 	rhs_ijk  = -1.0 * uarray[k][j][i];	
// 	f[k][j][i] = udot[k][j][i] - rhs_ijk;
//       }
//     }
//   }
//   /* Restore vectors */
//   DMDAVecRestoreArrayRead(da,localU,&uarray);
//   DMDAVecRestoreArray(da,Udot,&udot);
//   DMDAVecRestoreArray(da,F,&f);
  
//   // DMLocalToGlobalBegin( da , localF , INSERT_VALUES , F );
//   // DMLocalToGlobalEnd(   da , localF , INSERT_VALUES , F );
  
//   DMRestoreLocalVector(da,&localU);
//   // DMRestoreLocalVector(da,&localUdot);
//   // DMRestoreLocalVector(da,&localF);
  
//   PetscFunctionReturn(0);

// }


// PetscErrorCode FormRHSTEST(TS ts,PetscReal t,Vec U,Vec F,void *ctx) {

//   // Computes F = RHSfunction
  
//   PetscErrorCode ierr;
//   AppCtx         *user=(AppCtx*)ctx;
//   DM             da   = (DM)user->da;
//   PetscInt       i,j,k,Mx,My,Mz,xs,ys,zs,xm,ym,zm;
//   PetscScalar      hx,hy,hz,sx,sy,sz;
//   PetscScalar    u,***uarray,***f,***eps_2_array, ***sigma_array;
//   Vec            localU, localF, local_eps_2, local_sigma;
//   PetscScalar l_i,l_ip1,l_im1,l_jp1,l_jm1,l_km1,l_kp1;
//   PetscScalar dxx,dyy,dzz;
//   PetscScalar q_im1,q_ip1,q_jm1,q_jp1,q_km1,q_kp1,q_0;
//   PetscScalar sigma,m,eps_2;

//   PetscFunctionBeginUser;
  
//   DMGetLocalVector(da,&localU);
//   DMGetLocalVector(da,&localF);
  
//   DMDAGetInfo( da ,
//                PETSC_IGNORE,
//                &Mx , &My , &Mz ,
//                PETSC_IGNORE , PETSC_IGNORE , PETSC_IGNORE ,
//                PETSC_IGNORE ,
//                PETSC_IGNORE ,
//                PETSC_IGNORE , PETSC_IGNORE , PETSC_IGNORE ,
//                PETSC_IGNORE );
  
//   DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU);
//   DMGlobalToLocalEnd(da,U,INSERT_VALUES,localU);
//   DMGlobalToLocalBegin(da,F,INSERT_VALUES,localF);
//   DMGlobalToLocalEnd(da,F,INSERT_VALUES,localF);
  
//   DMDAVecGetArrayRead(da,localU,&uarray);
//   DMDAVecGetArray(da,localF,&f);  
  
//   DMDAGetCorners( da ,
//                   &xs , &ys , &zs ,
//                   &xm , &ym , &zm );
  
//   /* Compute function over the locally owned part of the grid */
//   for ( k = zs ; k < zs+zm ; k++ ) {
//     for ( j = ys ; j < ys+ym ; j++ ) {
//       for ( i = xs ; i < xs+xm ; i++ ) {
// 	f[k][j][i]  = -1.0 * uarray[k][j][i];	
//       }
//     }
//   }

//   /* Restore vectors */
//   DMDAVecRestoreArrayRead(da,localU,&uarray);
//   DMDAVecRestoreArray(da,localF,&f);
  
//   DMLocalToGlobalBegin( da , localF , INSERT_VALUES , F );
//   DMLocalToGlobalEnd(   da , localF , INSERT_VALUES , F );
  
//   DMRestoreLocalVector(da,&localU);
//   DMRestoreLocalVector(da,&localF);
  
//   PetscFunctionReturn(0);

// }

// FUNCTIONS FOR THE SPLIT-CH SOLVER

PetscScalar*** FormLocalRHS_CH_split_c( DMDALocalInfo *info ,
                                        PetscScalar ***uarray ,
                                        PetscScalar ***phiarray ,
                                        PetscScalar ***rhs ,
                                        PetscScalar ***sigma_array ,
                                        AppCtx *user ) {

  // Function to evaluate RHS of CH dynamics for a given local process
  // Output: rhs array, containing RHS evaluations on local process
  
  // NOTE: these CH eqns are dimensionless with domain length scale = 1. Physical domain size shows up in L_omega.
  PetscScalar hx = 1.0 / (PetscReal)(info->mx-1);
  PetscScalar sx = 1.0/(hx*hx);
  PetscScalar hy = 1.0 / (PetscReal)(info->my-1);
  PetscScalar sy = 1.0/(hy*hy);
  PetscScalar hz = 1.0 / (PetscReal)(info->mz-1);
  PetscScalar sz = 1.0/(hz*hz);

  // Set boundary node function depending on user-defined BCs
  void (*set_boundary_ghost_nodes)( AppCtx* , PetscScalar*** , PetscInt , PetscInt , PetscInt , PetscInt , PetscInt , PetscInt );
  if (user->boundary.compare("dirichlet") == 0) {
    set_boundary_ghost_nodes = set_boundary_ghost_nodes_dirichlet_singleframe;
  }
  else {
    // NEED TO IMPLEMENT THIS
  }

  /* Compute function over the locally owned part of the grid */
  for (int k=info->zs; k<info->zs+info->zm; k++) {
    for (int j=info->ys; j<info->ys+info->ym; j++) {
      for (int i=info->xs; i<info->xs+info->xm; i++) {

        (*set_boundary_ghost_nodes)( user , uarray , info->mx , info->my , info->mz , i , j , k );
      
        // dc/dt = laplacian( phi ) - sigma*(c - m)
        PetscScalar sigma = sigma_array[k][j][i];
        PetscScalar m     = user->m;
      
        // Term: laplacian( phi )
        PetscScalar l_i     = phiarray[k][j][i];
        PetscScalar l_im1   = phiarray[k][j][i-1];
        PetscScalar l_ip1   = phiarray[k][j][i+1];
        PetscScalar l_jm1   = phiarray[k][j-1][i];
        PetscScalar l_jp1   = phiarray[k][j+1][i];
        PetscScalar l_km1   = phiarray[k-1][j][i];
        PetscScalar l_kp1   = phiarray[k+1][j][i];

        PetscScalar dxx     = sx * ( l_ip1 + l_im1 - 2.0 * l_i );
        PetscScalar dyy     = sy * ( l_jp1 + l_jm1 - 2.0 * l_i );
        PetscScalar dzz     = sz * ( l_kp1 + l_km1 - 2.0 * l_i );

        rhs[k][j][i]  = dxx + dyy + dzz;
      
        // Term: -sigma*(c - m)
        rhs[k][j][i] += -sigma * ( uarray[k][j][i] - m );

      }
    }
  }
  
  return rhs;

}

PetscScalar*** FormLocalRHS_CH_split_phi( DMDALocalInfo *info ,
                                          PetscScalar ***uarray ,
                                          PetscScalar ***rhs ,
                                          PetscScalar ***eps2_array ,
                                          AppCtx *user ) {

  // Function to evaluate RHS of CH dynamics for a given local process
  // Output: rhs array, containing RHS evaluations on local process
  
  // NOTE: these CH eqns are dimensionless with domain length scale = 1. Physical domain size shows up in L_omega.
  PetscScalar hx = 1.0 / (PetscReal)(info->mx-1);
  PetscScalar sx = 1.0/(hx*hx);
  PetscScalar hy = 1.0 / (PetscReal)(info->my-1);
  PetscScalar sy = 1.0/(hy*hy);
  PetscScalar hz = 1.0 / (PetscReal)(info->mz-1);
  PetscScalar sz = 1.0/(hz*hz);

  // Set boundary node function depending on user-defined BCs
  void (*set_boundary_ghost_nodes)( AppCtx* , PetscScalar*** , PetscInt , PetscInt , PetscInt , PetscInt , PetscInt , PetscInt );
  if (user->boundary.compare("dirichlet") == 0) {
    set_boundary_ghost_nodes = set_boundary_ghost_nodes_dirichlet_singleframe;
  }
  else {
    // NEED TO IMPLEMENT THIS
  }

  /* Compute function over the locally owned part of the grid */
  for (int k=info->zs; k<info->zs+info->zm; k++) {
    for (int j=info->ys; j<info->ys+info->ym; j++) {
      for (int i=info->xs; i<info->xs+info->xm; i++) {

        (*set_boundary_ghost_nodes)( user , uarray , info->mx , info->my , info->mz , i , j , k );
      
        // dphi/dt = -eps2 * laplacian( u ) + ( u^3 - u )
      
        // Term: -eps2 * laplacian( u )
        PetscScalar l_i     = uarray[k][j][i];
        PetscScalar l_im1   = uarray[k][j][i-1];
        PetscScalar l_ip1   = uarray[k][j][i+1];
        PetscScalar l_jm1   = uarray[k][j-1][i];
        PetscScalar l_jp1   = uarray[k][j+1][i];
        PetscScalar l_km1   = uarray[k-1][j][i];
        PetscScalar l_kp1   = uarray[k+1][j][i];
      
        PetscScalar dxx     = sx * ( l_ip1 + l_im1 - 2.0 * l_i );
        PetscScalar dyy     = sy * ( l_jp1 + l_jm1 - 2.0 * l_i );
        PetscScalar dzz     = sz * ( l_kp1 + l_km1 - 2.0 * l_i );

        rhs[k][j][i]        = -eps2_array[k][j][i] * ( dxx + dyy + dzz );
      
        // Term: ( u^3 - u )
        rhs[k][j][i]       += uarray[k][j][i] * uarray[k][j][i] * uarray[k][j][i] - uarray[k][j][i];
      
      }
    }
  }
  
  return rhs;

}

PetscErrorCode FormIFunction_CH_split(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx) {

  // Computes residual F = Udot - RHSFunction

  AppCtx         *user = (AppCtx*)ctx;
  DM             pack  = user->pack;
  DM             da_c , da_phi;
  DMDALocalInfo  info_c , info_phi;
  PetscScalar    ***carray,***phiarray, ***f_c,***f_phi, ***udot_c,***udot_phi, ***eps_2_array,***sigma_array, ***rhs_c, ***rhs_phi;
  Vec            local_c,local_phi, local_eps_2,local_sigma, local_crhs,local_phirhs, U_c,U_phi, Udot_c,Udot_phi, F_c , F_phi;
  Vec            local_udotC,local_udotPhi, local_fC,local_fPhi;
  
  PetscFunctionBeginUser;

  DMCompositeGetEntries( pack , &da_c    , &da_phi );
  DMCompositeGetAccess(  pack , U        , &U_c      , &U_phi );
  DMCompositeGetAccess(  pack , Udot     , &Udot_c   , &Udot_phi );
  DMCompositeGetAccess(  pack , F        , &F_c      , &F_phi );
  
  DMGetLocalVector( da_c   , &local_c);
  DMGetLocalVector( da_phi , &local_phi);
  DMGetLocalVector( da_c   , &local_crhs);
  DMGetLocalVector( da_phi , &local_phirhs);
  DMGetLocalVector( da_phi , &local_eps_2);
  DMGetLocalVector( da_c   , &local_sigma);
  DMGetLocalVector( da_c   , &local_udotC);
  DMGetLocalVector( da_c   , &local_fC);
  DMGetLocalVector( da_phi , &local_udotPhi);
  DMGetLocalVector( da_phi , &local_fPhi);

  DMDAGetLocalInfo( da_c , &info_c );
  DMDAGetLocalInfo( da_c , &info_phi );
  
  DMGlobalToLocalBegin( da_c , U_c , INSERT_VALUES , local_c );
  DMGlobalToLocalEnd(   da_c , U_c , INSERT_VALUES , local_c );
  DMGlobalToLocalBegin( da_c , Udot_c , INSERT_VALUES , local_udotC );
  DMGlobalToLocalEnd(   da_c , Udot_c , INSERT_VALUES , local_udotC );
  DMGlobalToLocalBegin( da_c , F_c , INSERT_VALUES , local_fC );
  DMGlobalToLocalEnd(   da_c , F_c , INSERT_VALUES , local_fC );
  DMGlobalToLocalBegin( da_c , user->sigma , INSERT_VALUES , local_sigma );
  DMGlobalToLocalEnd(   da_c , user->sigma , INSERT_VALUES , local_sigma );
  DMGlobalToLocalBegin( da_phi , U_phi , INSERT_VALUES , local_phi );
  DMGlobalToLocalEnd(   da_phi , U_phi , INSERT_VALUES , local_phi );
  DMGlobalToLocalBegin( da_phi , Udot_phi , INSERT_VALUES , local_udotPhi );
  DMGlobalToLocalEnd(   da_phi , Udot_phi , INSERT_VALUES , local_udotPhi );
  DMGlobalToLocalBegin( da_phi , F_phi , INSERT_VALUES , local_fPhi );
  DMGlobalToLocalEnd(   da_phi , F_phi , INSERT_VALUES , local_fPhi );
  DMGlobalToLocalBegin( da_phi , user->eps_2 , INSERT_VALUES , local_eps_2 );
  DMGlobalToLocalEnd(   da_phi , user->eps_2 , INSERT_VALUES , local_eps_2 );

  DMDAVecGetArray(     da_c , local_c , &carray );
  DMDAVecGetArrayRead( da_c , local_sigma , &sigma_array );
  DMDAVecGetArrayRead( da_c , local_udotC , &udot_c );
  DMDAVecGetArray(     da_c , local_fC , &f_c );
  DMDAVecGetArray(     da_c , local_crhs , &rhs_c );
  DMDAVecGetArray(     da_phi , local_phi   , &phiarray );
  DMDAVecGetArrayRead( da_phi , local_eps_2 , &eps_2_array );
  DMDAVecGetArrayRead( da_phi , local_udotPhi , &udot_phi );
  DMDAVecGetArray(     da_phi , local_fPhi , &f_phi );
  DMDAVecGetArray(     da_phi , local_phirhs , &rhs_phi );
  
  /* Compute function over the locally owned part of the grid */
  rhs_c       = FormLocalRHS_CH_split_c(    &info_c   , carray   , phiarray , rhs_c , sigma_array , user );
  f_c         = FormLocalResidual(          &info_c   , carray   , f_c      , udot_c , rhs_c , user );
  rhs_phi     = FormLocalRHS_CH_split_phi(  &info_phi , carray   , rhs_phi  , eps_2_array , user );
  f_phi       = FormLocalResidual(          &info_phi , phiarray , f_phi    , phiarray    , rhs_phi , user );
  
  /* Restore vectors */
  DMDAVecRestoreArray(da_c,local_c,&carray);
  DMDAVecRestoreArrayRead(da_c,local_sigma,&sigma_array);
  DMDAVecRestoreArrayRead(da_c,local_udotC,&udot_c);
  DMDAVecRestoreArray(da_c,local_crhs,&rhs_c);
  DMDAVecRestoreArrayRead(da_phi,local_eps_2,&eps_2_array);
  DMDAVecRestoreArray(da_phi,local_phi,&phiarray);
  DMDAVecRestoreArrayRead(da_phi,local_udotPhi,&udot_phi);
  DMDAVecRestoreArray(da_phi,local_phirhs,&rhs_phi);

  DMDAVecRestoreArray(da_c  ,local_fC  ,&f_c);
  DMDAVecRestoreArray(da_phi,local_fPhi,&f_phi);

  DMLocalToGlobalBegin( da_c   , local_fC   , INSERT_VALUES , F_c );
  DMLocalToGlobalEnd(   da_c   , local_fC   , INSERT_VALUES , F_c );
  DMLocalToGlobalBegin( da_phi , local_fPhi , INSERT_VALUES , F_phi );
  DMLocalToGlobalEnd(   da_phi , local_fPhi , INSERT_VALUES , F_phi );

  DMCompositeRestoreAccess(  pack , U        , &U_c    , &U_phi );
  DMCompositeRestoreAccess(  pack , Udot     , &Udot_c , &Udot_phi );
  DMCompositeRestoreAccess(  pack , F        , &F_c    , &F_phi );
  
  DMRestoreLocalVector(da_c,&local_c);
  DMRestoreLocalVector(da_c,&local_crhs);
  DMRestoreLocalVector(da_c,&local_sigma);
  DMRestoreLocalVector(da_phi,&local_phi);
  DMRestoreLocalVector(da_phi,&local_phirhs);
  DMRestoreLocalVector(da_phi,&local_eps_2);

  DMRestoreLocalVector(da_c,&local_udotC);
  DMRestoreLocalVector(da_c,&local_fC);
  DMRestoreLocalVector(da_phi,&local_udotPhi);
  DMRestoreLocalVector(da_phi,&local_fPhi);
  
  PetscFunctionReturn(0);
  
}
