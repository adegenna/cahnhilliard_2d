#include <petscdmcomposite.h>
#include "rhs_implicit.h"
#include "boundary_conditions.h"
#include "temperature_dependence.h"

PetscErrorCode FormLocal_CH( DMDALocalInfo *info ,
                             PetscScalar **uarray ,
                             PetscScalar **eps_2_array ,
                             PetscScalar **sigma_array ,
                             PetscScalar **f , 
                             PetscScalar** udot ,
                             AppCtx *user ) {
    
  // NOTE: these CH eqns are dimensionless with domain length scale = 1. Physical domain size shows up in L_omega.
  PetscScalar hx = 1.0 / (PetscReal)(info->mx-1);
  PetscScalar sx = 1.0/(hx*hx);
  PetscScalar hy = 1.0 / (PetscReal)(info->my-1);
  PetscScalar sy = 1.0/(hy*hy);

  /* Compute function over the locally owned part of the grid */
  for (int j=info->ys; j<info->ys+info->ym; j++) {
    for (int i=info->xs; i<info->xs+info->xm; i++) {

      set_boundary_ghost_nodes( user , uarray , info->mx , info->my , i , j );
      
      /* Boundary conditions */
      ThirteenPointStencil stencil      = get_thirteen_point_stencil( user , uarray      , info->mx , info->my , i , j );
      ThirteenPointStencil stencil_eps2 = get_thirteen_point_stencil( user , eps_2_array , info->mx , info->my , i , j );
      
      // dc/dt = laplacian( c^3 - c ) - eps_2*biharm(c) - sigma*(c - m)
      PetscScalar sigma = sigma_array[j][i];
      PetscScalar m     = user->m;
      
      // Term: laplacian( c^3 - c )
      PetscScalar l_i     = stencil.c_i*stencil.c_i*stencil.c_i       - stencil.c_i;
      PetscScalar l_im1   = stencil.c_im1*stencil.c_im1*stencil.c_im1 - stencil.c_im1;
      PetscScalar l_ip1   = stencil.c_ip1*stencil.c_ip1*stencil.c_ip1 - stencil.c_ip1;
      PetscScalar l_jm1   = stencil.c_jm1*stencil.c_jm1*stencil.c_jm1 - stencil.c_jm1;
      PetscScalar l_jp1   = stencil.c_jp1*stencil.c_jp1*stencil.c_jp1 - stencil.c_jp1;

      PetscScalar dxx     = sx * ( l_ip1 + l_im1 - 2.0 * l_i );
      PetscScalar dyy     = sy * ( l_jp1 + l_jm1 - 2.0 * l_i );
      
      PetscScalar rhs_ij  = dxx + dyy;

      // Term: laplacian( -eps_2 * laplacian( c ) ) = laplacian( q )
      PetscScalar q_im1   = ( -stencil_eps2.c_im1 ) * ( sx * ( stencil.c_im2 + stencil.c_i   - 2.0*stencil.c_im1 ) + sy * ( stencil.c_ul  + stencil.c_bl  - 2.0*stencil.c_im1 ) );
      PetscScalar q_ip1   = ( -stencil_eps2.c_ip1 ) * ( sx * ( stencil.c_i   + stencil.c_ip2 - 2.0*stencil.c_ip1 ) + sy * ( stencil.c_ur  + stencil.c_br  - 2.0*stencil.c_ip1 ) );
      PetscScalar q_jm1   = ( -stencil_eps2.c_jm1 ) * ( sx * ( stencil.c_ul  + stencil.c_ur  - 2.0*stencil.c_jm1 ) + sy * ( stencil.c_jm2 + stencil.c_i   - 2.0*stencil.c_jm1 ) );
      PetscScalar q_jp1   = ( -stencil_eps2.c_jp1 ) * ( sx * ( stencil.c_bl  + stencil.c_br  - 2.0*stencil.c_jp1 ) + sy * ( stencil.c_jp2 + stencil.c_i   - 2.0*stencil.c_jp1 ) );
      PetscScalar q_0     = ( -stencil_eps2.c_i   ) * ( sx * ( stencil.c_im1 + stencil.c_ip1 - 2.0*stencil.c_i   ) + sy * ( stencil.c_jm1 + stencil.c_jp1 - 2.0*stencil.c_i   ) );
      
      rhs_ij += sx * ( q_im1 + q_ip1 - 2.0*q_0 ) + sy * ( q_jm1 + q_jp1 - 2.0*q_0 ); // laplacian( q )
      
      // Term: -sigma*(c - m)
      rhs_ij += -sigma * ( stencil.c_i - m );

      // Form f
      if ( user->boundary == 1 ) // Neumann: reset residuals explicitly 
        f[j][i] = reset_boundary_residual_values_for_neumann_bc( uarray , rhs_ij , udot[j][i] , info->mx , info->my , i , j );

      else if ( user->boundary == 3 ) // Bottom dirichlet, rest Neumann
        f[j][i] = reset_boundary_residual_values_for_dirichlet_bottom_neumann_remainder_bc( uarray , rhs_ij , udot[j][i] , info->mx , info->my , i , j );

      else if ( user->boundary == 4 ) // Bottom/top dirichlet, rest Neumann
        f[j][i] = reset_boundary_residual_values_for_dirichlet_topandbottom_neumann_remainder_bc( uarray , rhs_ij , udot[j][i] , info->mx , info->my , i , j );
      
      else // Dirichlet or periodic: just compute with ghost nodes
        f[j][i] = udot[j][i] - rhs_ij;

    }

  }

  PetscFunctionReturn(0);

}

PetscScalar** FormLocalRHS_thermal( DMDALocalInfo *info ,
                                    PetscScalar **Tarray ,
                                    PetscScalar **rhs ,
                                    PetscScalar **Tsource ,
                                    AppCtx *user ) {

  // Function to evaluate RHS of thermal dynamics for a given local process
  // Output: rhs array, containing RHS evaluations on local process

  // NOTE: these CH eqns are dimensionless with domain length scale = 1. Physical domain size shows up in L_omega.
  PetscScalar hx = 1.0 / (PetscReal)(info->mx-1);
  PetscScalar sx = 1.0 / (hx*hx);
  PetscScalar hy = 1.0 / (PetscReal)(info->my-1);
  PetscScalar sy = 1.0 / (hy*hy);
  
  /* Compute function over the locally owned part of the grid */
  for (int j=info->ys; j<info->ys+info->ym; j++) {
    for (int i=info->xs; i<info->xs+info->xm; i++) {

      set_boundary_ghost_nodes( user , Tarray  , info->mx , info->my , i , j );
      //set_boundary_ghost_nodes( user , Tsource , info->mx , info->my , i , j );

      // dT/dt = D_T * laplacian( T ) + S
      
      // Term: D_T * laplacian( T )
      PetscScalar dxx     = sx * ( Tarray[j][i+1] + Tarray[j][i-1] - 2.0 * Tarray[j][i] );
      PetscScalar dyy     = sy * ( Tarray[j+1][i] + Tarray[j-1][i] - 2.0 * Tarray[j][i] );
      
      rhs[j][i]           = user->D_T * ( dxx + dyy ) + Tsource[j][i];
      
    }

  }

  return rhs;
  
}

PetscScalar** FormLocal_thermal( DMDALocalInfo* info ,
                                 PetscScalar** Tarray ,
                                 PetscScalar** f , 
                                 PetscScalar** udot ,
                                 PetscScalar** Tsource ,
                                 PetscScalar** rhs ,
                                 AppCtx* user ) {
  
  // Apply correction to boundary fluxes and set residuals
  for (int j=info->ys; j<info->ys+info->ym; j++) {
    for (int i=info->xs; i<info->xs+info->xm; i++) {
      f[j][i] = reset_boundary_residual_values_for_neumann_bc( Tarray , rhs[j][i] , udot[j][i] , info->mx , info->my , i , j );
    }
  }

  return f;

}

PetscErrorCode FormIFunction_CH(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx) {

  // Computes residual F = Udot - RHSFunction

  PetscErrorCode ierr;
  AppCtx         *user = (AppCtx*)ctx;
  DM             da_c  = (DM)user->da_c;
  DMDALocalInfo  info_c;
  PetscScalar    u,**carray,**f,**udot, **eps_2_array, **sigma_array;
  Vec            local_c, local_eps_2, local_sigma, local_udot , local_f;
  
  PetscFunctionBeginUser;
  
  DMGetLocalVector( da_c , &local_c);
  DMGetLocalVector( da_c , &local_eps_2);
  DMGetLocalVector( da_c , &local_sigma);

  DMDAGetLocalInfo( da_c , &info_c );
  
  DMGlobalToLocalBegin( da_c , U , INSERT_VALUES , local_c );
  DMGlobalToLocalEnd(   da_c , U , INSERT_VALUES , local_c );
  DMGlobalToLocalBegin( da_c , Udot , INSERT_VALUES , local_udot );
  DMGlobalToLocalEnd(   da_c , Udot , INSERT_VALUES , local_udot );
  DMGlobalToLocalBegin( da_c , F , INSERT_VALUES , local_f );
  DMGlobalToLocalEnd(   da_c , F , INSERT_VALUES , local_f );
  DMGlobalToLocalBegin( da_c , user->eps_2 , INSERT_VALUES , local_eps_2 );
  DMGlobalToLocalEnd(   da_c , user->eps_2 , INSERT_VALUES , local_eps_2 );
  DMGlobalToLocalBegin( da_c , user->sigma , INSERT_VALUES , local_sigma );
  DMGlobalToLocalEnd(   da_c , user->sigma , INSERT_VALUES , local_sigma );

  DMDAVecGetArrayRead( da_c , local_c , &carray );
  DMDAVecGetArrayRead( da_c , local_eps_2 , &eps_2_array );
  DMDAVecGetArrayRead( da_c , local_sigma , &sigma_array );
  DMDAVecGetArrayRead( da_c , local_udot , &udot );
  DMDAVecGetArray(     da_c , local_f , &f );
  
  /* Compute function over the locally owned part of the grid */
  FormLocal_CH( &info_c , carray , eps_2_array , sigma_array , f , udot , user );

  /* Restore vectors */
  DMDAVecRestoreArrayRead(da_c,local_c,&carray);
  DMDAVecRestoreArrayRead(da_c,local_eps_2,&eps_2_array);
  DMDAVecRestoreArrayRead(da_c,local_sigma,&sigma_array);
  DMDAVecRestoreArrayRead(da_c,local_udot,&udot);
  DMDAVecRestoreArray(da_c,local_f,&f);

  DMLocalToGlobalBegin( da_c , local_f , INSERT_VALUES , F );
  DMLocalToGlobalEnd(   da_c , local_f , INSERT_VALUES , F );

  DMRestoreLocalVector(da_c,&local_c);
  DMRestoreLocalVector(da_c,&local_eps_2);
  DMRestoreLocalVector(da_c,&local_sigma);
  DMRestoreLocalVector(da_c,&local_udot);
  DMRestoreLocalVector(da_c,&local_f);
  
  
  PetscFunctionReturn(0);
  
}

PetscErrorCode FormRHS_thermal(TS ts,PetscReal t,Vec U,Vec F,void *ctx) {

  // Computes F = RHSfunction

  AppCtx         *user = (AppCtx*)ctx;
  DM             da_T  = (DM)user->da_T;
  DMDALocalInfo  info_T;
  PetscScalar    u,**Tarray,**Tsource,**rhs_thermal;
  Vec            local_T,local_Tsource,local_rhs;
  
  PetscFunctionBeginUser;
  
  DMGetLocalVector( da_T , &local_T );
  DMGetLocalVector( da_T , &local_Tsource );
  DMGetLocalVector( da_T , &local_rhs );
  
  DMDAGetLocalInfo( da_T , &info_T );
  
  DMGlobalToLocalBegin( da_T , U , INSERT_VALUES , local_T );
  DMGlobalToLocalEnd(   da_T , U , INSERT_VALUES , local_T );
  DMGlobalToLocalBegin( da_T , F , INSERT_VALUES , local_rhs );
  DMGlobalToLocalEnd(   da_T , F , INSERT_VALUES , local_rhs );
  DMGlobalToLocalBegin( da_T , user->temperature_source , INSERT_VALUES , local_Tsource );
  DMGlobalToLocalEnd(   da_T , user->temperature_source , INSERT_VALUES , local_Tsource );

  DMDAVecGetArrayRead( da_T , local_T , &Tarray );
  DMDAVecGetArrayRead( da_T , local_Tsource , &Tsource );
  DMDAVecGetArray(     da_T , local_rhs , &rhs_thermal );
  
  /* Compute function over the locally owned part of the grid */
  rhs_thermal = FormLocalRHS_thermal( &info_T , Tarray , rhs_thermal , Tsource , user );

  /* Restore vectors */
  DMDAVecRestoreArrayRead(da_T,local_T,&Tarray);
  DMDAVecRestoreArray(da_T,local_rhs,&rhs_thermal);

  DMLocalToGlobalBegin( da_T , local_rhs , INSERT_VALUES , F );
  DMLocalToGlobalEnd(   da_T , local_rhs , INSERT_VALUES , F );
  
  DMRestoreLocalVector(da_T,&local_rhs);
  DMRestoreLocalVector(da_T,&local_T);
  DMRestoreLocalVector(da_T,&local_Tsource);
  DMRestoreLocalVector(da_T,&local_rhs);

  PetscFunctionReturn(0);
  
}

PetscErrorCode FormIFunction_thermal(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx) {

  // Computes residual F = Udot - RHSFunction

  AppCtx         *user = (AppCtx*)ctx;
  DM             da_T  = (DM)user->da_T;
  DMDALocalInfo  info_T;
  PetscScalar    u,**Tarray,**Tsource,**f,**udot,**rhs_thermal;
  Vec            local_T,local_Tsource,local_Trhs , local_udot , local_f;
  
  PetscFunctionBeginUser;
  
  DMGetLocalVector( da_T , &local_T );
  DMGetLocalVector( da_T , &local_Tsource );
  DMGetLocalVector( da_T , &local_Trhs );
  DMGetLocalVector( da_T , &local_udot );
  DMGetLocalVector( da_T , &local_f );

  DMDAGetLocalInfo( da_T , &info_T );
  
  DMGlobalToLocalBegin( da_T , U , INSERT_VALUES , local_T );
  DMGlobalToLocalEnd(   da_T , U , INSERT_VALUES , local_T );
  DMGlobalToLocalBegin( da_T , Udot , INSERT_VALUES , local_udot );
  DMGlobalToLocalEnd(   da_T , Udot , INSERT_VALUES , local_udot );
  DMGlobalToLocalBegin( da_T , F , INSERT_VALUES , local_f );
  DMGlobalToLocalEnd(   da_T , F , INSERT_VALUES , local_f );
  DMGlobalToLocalBegin( da_T , user->temperature_source , INSERT_VALUES , local_Tsource );
  DMGlobalToLocalEnd(   da_T , user->temperature_source , INSERT_VALUES , local_Tsource );

  DMDAVecGetArrayRead( da_T , local_T , &Tarray );
  DMDAVecGetArrayRead( da_T , local_udot , &udot );
  DMDAVecGetArrayRead( da_T , local_Tsource , &Tsource );
  DMDAVecGetArray(     da_T , local_f , &f );
  DMDAVecGetArray(     da_T , local_Trhs , &rhs_thermal );
  
  /* Compute function over the locally owned part of the grid */
  rhs_thermal = FormLocalRHS_thermal( &info_T , Tarray , rhs_thermal , Tsource , user );
  f           = FormLocal_thermal(    &info_T , Tarray , f , udot , Tsource , rhs_thermal , user );

  /* Restore vectors */
  DMDAVecRestoreArrayRead(da_T,local_T,&Tarray);
  DMDAVecRestoreArrayRead(da_T,local_udot,&udot);  
  DMDAVecRestoreArray(da_T,local_f,&f);

  DMLocalToGlobalBegin( da_T , local_f , INSERT_VALUES , F );
  DMLocalToGlobalEnd(   da_T , local_f , INSERT_VALUES , F );
  
  DMRestoreLocalVector(da_T,&local_T);
  DMRestoreLocalVector(da_T,&local_Tsource);
  DMRestoreLocalVector(da_T,&local_Trhs);
  DMRestoreLocalVector(da_T,&local_f);
  DMRestoreLocalVector(da_T,&local_udot);
  
  PetscFunctionReturn(0);
  
}

PetscErrorCode FormIFunction_CH_coupled(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx) {

  // Computes residual F = Udot - RHSFunction

  PetscErrorCode ierr;
  AppCtx         *user = (AppCtx*)ctx;
  DM             pack  = user->pack;
  DM             da_c , da_T;
  DMDALocalInfo  info_c , info_T;
  PetscScalar    u,**carray,**Tarray,**Tsource,**f_c,**f_T,**udot, **udot_c, **udot_T, **eps_2_array, **sigma_array, **rhs_thermal;
  Vec            local_c, local_T, local_Tsource, local_eps_2, local_sigma, local_Trhs, U_c , U_T , Udot_c , Udot_T , F_c , F_T;
  Vec            local_udotC , local_udotT , local_fC , local_fT;
  
  PetscFunctionBeginUser;

  DMCompositeGetEntries( pack , &da_c    , &da_T );
  DMCompositeGetAccess(  pack , U        , &U_c , &U_T );
  DMCompositeGetAccess(  pack , Udot     , &Udot_c , &Udot_T );
  DMCompositeGetAccess(  pack , F , &F_c , &F_T );
  
  DMGetLocalVector( da_c , &local_c);
  DMGetLocalVector( da_T , &local_T);
  DMGetLocalVector( da_T , &local_Tsource);
  DMGetLocalVector( da_T , &local_Trhs);
  DMGetLocalVector( da_c , &local_eps_2);
  DMGetLocalVector( da_c , &local_sigma);
  DMGetLocalVector( da_c , &local_udotC);
  DMGetLocalVector( da_c , &local_fC);
  DMGetLocalVector( da_T , &local_udotT);
  DMGetLocalVector( da_T , &local_fT);  

  DMDAGetLocalInfo( da_c , &info_c );
  DMDAGetLocalInfo( da_T , &info_T );
  
  DMGlobalToLocalBegin( da_c , U_c , INSERT_VALUES , local_c );
  DMGlobalToLocalEnd(   da_c , U_c , INSERT_VALUES , local_c );
  DMGlobalToLocalBegin( da_c , user->eps_2 , INSERT_VALUES , local_eps_2 );
  DMGlobalToLocalEnd(   da_c , user->eps_2 , INSERT_VALUES , local_eps_2 );
  DMGlobalToLocalBegin( da_c , user->sigma , INSERT_VALUES , local_sigma );
  DMGlobalToLocalEnd(   da_c , user->sigma , INSERT_VALUES , local_sigma );
  DMGlobalToLocalBegin( da_T , U_T , INSERT_VALUES , local_T );
  DMGlobalToLocalEnd(   da_T , U_T , INSERT_VALUES , local_T );
  DMGlobalToLocalBegin( da_T , user->temperature_source , INSERT_VALUES , local_Tsource );
  DMGlobalToLocalEnd(   da_T , user->temperature_source , INSERT_VALUES , local_Tsource );
  DMGlobalToLocalBegin( da_c , Udot_c , INSERT_VALUES , local_udotC );
  DMGlobalToLocalEnd(   da_c , Udot_c , INSERT_VALUES , local_udotC );
  DMGlobalToLocalBegin( da_c , F_c , INSERT_VALUES , local_fC );
  DMGlobalToLocalEnd(   da_c , F_c , INSERT_VALUES , local_fC );
  DMGlobalToLocalBegin( da_T , Udot_T , INSERT_VALUES , local_udotT );
  DMGlobalToLocalEnd(   da_T , Udot_T , INSERT_VALUES , local_udotT );
  DMGlobalToLocalBegin( da_T , F_T , INSERT_VALUES , local_fT );
  DMGlobalToLocalEnd(   da_T , F_T , INSERT_VALUES , local_fT );

  DMDAVecGetArray( da_c , local_c , &carray );
  DMDAVecGetArrayRead( da_c , local_eps_2 , &eps_2_array );
  DMDAVecGetArrayRead( da_c , local_sigma , &sigma_array );
  DMDAVecGetArrayRead( da_c , local_udotC , &udot_c );
  DMDAVecGetArray(     da_c , local_fC , &f_c );
  DMDAVecGetArrayRead( da_T , local_udotT , &udot_T );
  DMDAVecGetArray(     da_T , local_T , &Tarray );
  DMDAVecGetArrayRead( da_T , local_Tsource , &Tsource );
  DMDAVecGetArray(     da_T , local_Trhs , &rhs_thermal );
  DMDAVecGetArray(     da_T , local_fT , &f_T );
  
  /* Compute function over the locally owned part of the grid */
  compute_eps2_and_sigma_from_temperature( user , U );
  FormLocal_CH(      &info_c , carray , eps_2_array , sigma_array , f_c , udot_c , user );
  rhs_thermal = FormLocalRHS_thermal( &info_T , Tarray , rhs_thermal , Tsource , user );
  f_T         = FormLocal_thermal(    &info_T , Tarray , f_T , udot_T , Tsource , rhs_thermal , user );
  
  /* Restore vectors */
  DMDAVecRestoreArray(da_c,local_c,&carray);
  DMDAVecRestoreArrayRead(da_c,local_eps_2,&eps_2_array);
  DMDAVecRestoreArrayRead(da_c,local_sigma,&sigma_array);
  DMDAVecRestoreArrayRead(da_c,local_udotC,&udot_c);
  DMDAVecRestoreArray(da_T,local_T,&Tarray);
  DMDAVecRestoreArrayRead(da_T,local_Tsource,&Tsource);
  DMDAVecRestoreArray(da_T,local_Trhs,&rhs_thermal);
  DMDAVecRestoreArrayRead(da_T,local_udotT,&udot_T);

  DMDAVecRestoreArray(da_c,local_fC,&f_c);
  DMDAVecRestoreArray(da_T,local_fT,&f_T);

  DMLocalToGlobalBegin( da_T , local_fT , INSERT_VALUES , F_T );
  DMLocalToGlobalEnd(   da_T , local_fT , INSERT_VALUES , F_T );
  DMLocalToGlobalBegin( da_c , local_fC , INSERT_VALUES , F_c );
  DMLocalToGlobalEnd(   da_c , local_fC , INSERT_VALUES , F_c );
  
  DMCompositeRestoreAccess(  pack , U        , &U_c , &U_T );
  DMCompositeRestoreAccess(  pack , Udot     , &Udot_c , &Udot_T );
  DMCompositeRestoreAccess(  pack , F        , &F_c , &F_T );
  
  DMRestoreLocalVector(da_c,&local_c);
  DMRestoreLocalVector(da_T,&local_T);
  DMRestoreLocalVector(da_T,&local_Tsource);
  DMRestoreLocalVector(da_T,&local_Trhs);
  DMRestoreLocalVector(da_c,&local_eps_2);
  DMRestoreLocalVector(da_c,&local_sigma);

  DMRestoreLocalVector(da_c,&local_udotC);
  DMRestoreLocalVector(da_c,&local_fC);
  DMRestoreLocalVector(da_T,&local_udotT);
  DMRestoreLocalVector(da_T,&local_fT);  
  
  PetscFunctionReturn(0);
  
}
