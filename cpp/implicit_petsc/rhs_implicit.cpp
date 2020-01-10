#include <petscdmcomposite.h>
#include "rhs_implicit.h"
#include "boundary_conditions.h"
#include "temperature_dependence.h"

PetscScalar** FormLocal_CH( DMDALocalInfo *info ,
                            PetscScalar **uarray ,
                            PetscScalar **f , 
                            PetscScalar **udot ,
                            PetscScalar **rhs ,
                            AppCtx *user ) {
  
  // Apply correction to boundary fluxes and set residuals
  for (int j=info->ys; j<info->ys+info->ym; j++) {
    for (int i=info->xs; i<info->xs+info->xm; i++) {

      f[j][i] = user->residualFunction( uarray , rhs[j][i] , udot[j][i] , info->mx , info->my , i , j );
      
    }
  }

  return f;

}

PetscScalar** FormLocalRHS_CH( DMDALocalInfo *info ,
                               PetscScalar **uarray ,
                               PetscScalar **rhs ,
                               PetscScalar **eps_2_array ,
                               PetscScalar **sigma_array ,
                               AppCtx *user ) {

  // Function to evaluate RHS of CH dynamics for a given local process
  // Output: rhs array, containing RHS evaluations on local process
  
  // NOTE: these CH eqns are dimensionless with domain length scale = 1. Physical domain size shows up in L_omega.
  PetscScalar hx = 1.0 / (PetscReal)(info->mx-1);
  PetscScalar sx = 1.0/(hx*hx);
  PetscScalar hy = 1.0 / (PetscReal)(info->my-1);
  PetscScalar sy = 1.0/(hy*hy);

  // Set boundary node function depending on user-defined BCs
  void (*set_boundary_ghost_nodes)( AppCtx* , PetscScalar** , PetscInt , PetscInt , PetscInt , PetscInt );
  if (user->boundary.compare("dirichlet") == 0) {
    set_boundary_ghost_nodes = set_boundary_ghost_nodes_dirichlet;
  }
  else if (user->boundary.compare("neumann") == 0) {
    set_boundary_ghost_nodes = set_boundary_ghost_nodes_neumann;
  }
  else {
    // For mixed neumann/dirichlet BCs, just set ghost nodes to dirichlet values and reset boundary fluxes later on
    // TODO: this does not work for periodic BCs (need to implement those)
    set_boundary_ghost_nodes = set_boundary_ghost_nodes_dirichlet;
  }

  /* Compute function over the locally owned part of the grid */
  for (int j=info->ys; j<info->ys+info->ym; j++) {
    for (int i=info->xs; i<info->xs+info->xm; i++) {

      (*set_boundary_ghost_nodes)( user , uarray , info->mx , info->my , i , j );
      
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
      
      rhs[j][i]  = dxx + dyy;

      // Term: laplacian( -eps_2 * laplacian( c ) ) = laplacian( q )
      PetscScalar q_im1   = ( -stencil_eps2.c_im1 ) * ( sx * ( stencil.c_im2 + stencil.c_i   - 2.0*stencil.c_im1 ) + sy * ( stencil.c_ul  + stencil.c_bl  - 2.0*stencil.c_im1 ) );
      PetscScalar q_ip1   = ( -stencil_eps2.c_ip1 ) * ( sx * ( stencil.c_i   + stencil.c_ip2 - 2.0*stencil.c_ip1 ) + sy * ( stencil.c_ur  + stencil.c_br  - 2.0*stencil.c_ip1 ) );
      PetscScalar q_jm1   = ( -stencil_eps2.c_jm1 ) * ( sx * ( stencil.c_ul  + stencil.c_ur  - 2.0*stencil.c_jm1 ) + sy * ( stencil.c_jm2 + stencil.c_i   - 2.0*stencil.c_jm1 ) );
      PetscScalar q_jp1   = ( -stencil_eps2.c_jp1 ) * ( sx * ( stencil.c_bl  + stencil.c_br  - 2.0*stencil.c_jp1 ) + sy * ( stencil.c_jp2 + stencil.c_i   - 2.0*stencil.c_jp1 ) );
      PetscScalar q_0     = ( -stencil_eps2.c_i   ) * ( sx * ( stencil.c_im1 + stencil.c_ip1 - 2.0*stencil.c_i   ) + sy * ( stencil.c_jm1 + stencil.c_jp1 - 2.0*stencil.c_i   ) );
      
      rhs[j][i] += sx * ( q_im1 + q_ip1 - 2.0*q_0 ) + sy * ( q_jm1 + q_jp1 - 2.0*q_0 ); // laplacian( q )
      
      // Term: -sigma*(c - m)
      rhs[j][i] += -sigma * ( stencil.c_i - m );

    }

  }

  return rhs;

}

PetscErrorCode FormRHS_CH(TS ts,PetscReal t,Vec U,Vec F,void *ctx) {

  // Computes F = RHSfunction

  AppCtx         *user = (AppCtx*)ctx;
  DM              pack = (DM)user->pack;
  DM              da_c , da_T;
  DMDALocalInfo   info_c;
  PetscScalar     **carray,**eps2,**sigma,**rhs_c;
  Vec             local_c,local_eps2,local_sigma,local_rhs;
  
  PetscFunctionBeginUser;

  DMCompositeGetEntries( pack , &da_c , &da_T );
  
  DMGetLocalVector( da_c , &local_c );
  DMGetLocalVector( da_c , &local_eps2 );
  DMGetLocalVector( da_c , &local_sigma );
  DMGetLocalVector( da_c , &local_rhs );
  
  DMDAGetLocalInfo( da_c , &info_c );
  
  DMGlobalToLocalBegin( da_c , U , INSERT_VALUES , local_c );
  DMGlobalToLocalEnd(   da_c , U , INSERT_VALUES , local_c );
  DMGlobalToLocalBegin( da_c , F , INSERT_VALUES , local_rhs );
  DMGlobalToLocalEnd(   da_c , F , INSERT_VALUES , local_rhs );
  DMGlobalToLocalBegin( da_c , user->eps_2 , INSERT_VALUES , local_eps2 );
  DMGlobalToLocalEnd(   da_c , user->eps_2 , INSERT_VALUES , local_eps2 );
  DMGlobalToLocalBegin( da_c , user->sigma , INSERT_VALUES , local_sigma );
  DMGlobalToLocalEnd(   da_c , user->sigma , INSERT_VALUES , local_sigma );
  
  DMDAVecGetArrayRead( da_c , local_c , &carray );
  DMDAVecGetArrayRead( da_c , local_eps2 , &eps2 );
  DMDAVecGetArrayRead( da_c , local_sigma , &sigma );
  DMDAVecGetArray(     da_c , local_rhs , &rhs_c );
  
  /* Compute function over the locally owned part of the grid */
  rhs_c = FormLocalRHS_CH( &info_c , carray , rhs_c , eps2 , sigma , user );

  /* Restore vectors */
  DMDAVecRestoreArrayRead( da_c , local_c , &carray );
  DMDAVecRestoreArrayRead( da_c , local_eps2 , &eps2 );
  DMDAVecRestoreArrayRead( da_c , local_sigma , &sigma );
  DMDAVecRestoreArray(     da_c , local_rhs , &rhs_c );

  DMLocalToGlobalBegin( da_c , local_rhs , INSERT_VALUES , F );
  DMLocalToGlobalEnd(   da_c , local_rhs , INSERT_VALUES , F );
  
  DMRestoreLocalVector( da_c , &local_c );
  DMRestoreLocalVector( da_c , &local_eps2 );
  DMRestoreLocalVector( da_c , &local_sigma );
  DMRestoreLocalVector( da_c , &local_rhs );

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

  // Set boundary node function depending on user-defined BCs
  void (*set_boundary_ghost_nodes)( AppCtx* , PetscScalar** , PetscInt , PetscInt , PetscInt , PetscInt );
  if (user->boundary.compare("dirichlet") == 0) {
    set_boundary_ghost_nodes = set_boundary_ghost_nodes_dirichlet_singleframe_thermal;
  }
  else if (user->boundary.compare("neumann") == 0) {
    set_boundary_ghost_nodes = set_boundary_ghost_nodes_neumann;
  }
  else {
    // For mixed neumann/dirichlet BCs, just set ghost nodes to dirichlet values and reset boundary fluxes later on
    // TODO: this does not work for periodic BCs (need to implement those)
    set_boundary_ghost_nodes = set_boundary_ghost_nodes_dirichlet;
  }
  
  /* Compute function over the locally owned part of the grid */
  for (int j=info->ys; j<info->ys+info->ym; j++) {
    for (int i=info->xs; i<info->xs+info->xm; i++) {

      (*set_boundary_ghost_nodes)( user , Tarray  , info->mx , info->my , i , j );

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

PetscErrorCode FormRHS_CH_coupled(TS ts,PetscReal t,Vec U,Vec F,void *ctx) {

  // Computes F = RHSfunction

  AppCtx         *user = (AppCtx*)ctx;
  DM              pack = (DM)user->pack;
  DM              da_c , da_T;
  DMDALocalInfo   info_c;
  PetscScalar     **carray,**eps2,**sigma,**rhs_c;
  Vec             local_c,local_eps2,local_sigma,local_rhs_c;
  DMDALocalInfo   info_T;
  PetscScalar     **Tarray,**Tsource,**rhs_thermal;
  Vec             local_T,local_Tsource,local_rhs_T;
  Vec             U_c , U_T , F_c , F_T;
  
  PetscFunctionBeginUser;

  // Get composite stuff
  DMCompositeGetEntries( pack , &da_c , &da_T );
  DMCompositeGetAccess(  pack , U     , &U_c , &U_T );
  DMCompositeGetAccess(  pack , F     , &F_c , &F_T );

  // Get CH data types
  DMGetLocalVector( da_c , &local_c );
  DMGetLocalVector( da_c , &local_eps2 );
  DMGetLocalVector( da_c , &local_sigma );
  DMGetLocalVector( da_c , &local_rhs_c );
  
  DMDAGetLocalInfo( da_c , &info_c );
  
  DMGlobalToLocalBegin( da_c , U_c , INSERT_VALUES , local_c );
  DMGlobalToLocalEnd(   da_c , U_c , INSERT_VALUES , local_c );
  DMGlobalToLocalBegin( da_c , F_c , INSERT_VALUES , local_rhs_c );
  DMGlobalToLocalEnd(   da_c , F_c , INSERT_VALUES , local_rhs_c );
  DMGlobalToLocalBegin( da_c , user->eps_2 , INSERT_VALUES , local_eps2 );
  DMGlobalToLocalEnd(   da_c , user->eps_2 , INSERT_VALUES , local_eps2 );
  DMGlobalToLocalBegin( da_c , user->sigma , INSERT_VALUES , local_sigma );
  DMGlobalToLocalEnd(   da_c , user->sigma , INSERT_VALUES , local_sigma );
  
  DMDAVecGetArrayRead( da_c , local_c , &carray );
  DMDAVecGetArrayRead( da_c , local_eps2 , &eps2 );
  DMDAVecGetArrayRead( da_c , local_sigma , &sigma );
  DMDAVecGetArray(     da_c , local_rhs_c , &rhs_c );

  // Get thermal data types
  DMGetLocalVector( da_T , &local_T );
  DMGetLocalVector( da_T , &local_Tsource );
  DMGetLocalVector( da_T , &local_rhs_T );
  
  DMDAGetLocalInfo( da_T , &info_T );
  
  DMGlobalToLocalBegin( da_T , U_T , INSERT_VALUES , local_T );
  DMGlobalToLocalEnd(   da_T , U_T , INSERT_VALUES , local_T );
  DMGlobalToLocalBegin( da_T , F_T , INSERT_VALUES , local_rhs_T );
  DMGlobalToLocalEnd(   da_T , F_T , INSERT_VALUES , local_rhs_T );
  DMGlobalToLocalBegin( da_T , user->temperature_source , INSERT_VALUES , local_Tsource );
  DMGlobalToLocalEnd(   da_T , user->temperature_source , INSERT_VALUES , local_Tsource );

  DMDAVecGetArrayRead( da_T , local_T , &Tarray );
  DMDAVecGetArrayRead( da_T , local_Tsource , &Tsource );
  DMDAVecGetArray(     da_T , local_rhs_T , &rhs_thermal );
  
  /* Compute function over the locally owned part of the grid */
  rhs_c       = FormLocalRHS_CH(      &info_c , carray , rhs_c , eps2 , sigma , user );
  rhs_thermal = FormLocalRHS_thermal( &info_T , Tarray , rhs_thermal , Tsource , user );
  
  /* Restore vectors */
  DMDAVecRestoreArrayRead( da_c , local_c , &carray );
  DMDAVecRestoreArrayRead( da_c , local_eps2 , &eps2 );
  DMDAVecRestoreArrayRead( da_c , local_sigma , &sigma );
  DMDAVecRestoreArray(     da_c , local_rhs_c , &rhs_c );
  DMDAVecRestoreArrayRead( da_T , local_T , &Tarray );
  DMDAVecRestoreArrayRead( da_T , local_Tsource , &Tsource );
  DMDAVecRestoreArray(     da_T , local_rhs_T , &rhs_thermal );
  
  DMLocalToGlobalBegin( da_c , local_rhs_c , INSERT_VALUES , F_c );
  DMLocalToGlobalEnd(   da_c , local_rhs_c , INSERT_VALUES , F_c );
  DMLocalToGlobalBegin( da_T , local_rhs_T , INSERT_VALUES , F_T );
  DMLocalToGlobalEnd(   da_T , local_rhs_T , INSERT_VALUES , F_T );

  DMCompositeRestoreAccess( pack , F , &F_c , &F_T );
  DMCompositeRestoreAccess( pack , U , &U_c , &U_T );

  DMRestoreLocalVector( da_c , &local_c );
  DMRestoreLocalVector( da_c , &local_eps2 );
  DMRestoreLocalVector( da_c , &local_sigma );
  DMRestoreLocalVector( da_c , &local_rhs_c );
  DMRestoreLocalVector( da_T , &local_T );
  DMRestoreLocalVector( da_T , &local_Tsource );
  DMRestoreLocalVector( da_T , &local_rhs_T );
  
  PetscFunctionReturn(0);
  
}

PetscErrorCode FormIFunction_CH(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx) {

  // Computes residual F = Udot - RHSFunction

  AppCtx         *user = (AppCtx*)ctx;
  DM             da_c  = (DM)user->da_c;
  DMDALocalInfo  info_c;
  PetscScalar    **carray,**f,**udot, **eps_2_array, **sigma_array, **rhs;
  Vec            local_c, local_eps_2, local_sigma, local_udot , local_f , local_rhs;
  
  PetscFunctionBeginUser;
  
  // Get CH data types
  DMGetLocalVector( da_c , &local_c );
  DMGetLocalVector( da_c , &local_eps_2 );
  DMGetLocalVector( da_c , &local_sigma );
  DMGetLocalVector( da_c , &local_rhs );
  DMGetLocalVector( da_c , &local_udot );
  DMGetLocalVector( da_c , &local_f );
  
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
  DMDAVecGetArray(     da_c , local_rhs , &rhs );

  /* Compute function over the locally owned part of the grid */
  rhs = FormLocalRHS_CH( &info_c , carray , rhs , eps_2_array , sigma_array , user );
  f   = FormLocal_CH(    &info_c , carray , f , udot , rhs , user );
  
  /* Restore vectors */
  DMDAVecRestoreArrayRead( da_c , local_c , &carray );
  DMDAVecRestoreArrayRead( da_c , local_eps_2 , &eps_2_array );
  DMDAVecRestoreArrayRead( da_c , local_sigma , &sigma_array );
  DMDAVecRestoreArrayRead( da_c , local_udot , &udot );
  DMDAVecRestoreArray(     da_c , local_f , &f );
  DMDAVecRestoreArray(     da_c , local_rhs , &rhs );

  DMLocalToGlobalBegin( da_c , local_f , INSERT_VALUES , F );
  DMLocalToGlobalEnd(   da_c , local_f , INSERT_VALUES , F );

  DMRestoreLocalVector( da_c , &local_c );
  DMRestoreLocalVector( da_c , &local_eps_2 );
  DMRestoreLocalVector( da_c , &local_sigma );
  DMRestoreLocalVector( da_c , &local_udot );
  DMRestoreLocalVector( da_c , &local_f );
  DMRestoreLocalVector( da_c , &local_rhs );
  
  PetscFunctionReturn(0);
  
}

PetscErrorCode FormRHS_thermal(TS ts,PetscReal t,Vec U,Vec F,void *ctx) {

  // Computes F = RHSfunction

  AppCtx         *user = (AppCtx*)ctx;
  DM             da_T  = (DM)user->da_T;
  DMDALocalInfo  info_T;
  PetscScalar    **Tarray,**Tsource,**rhs_thermal;
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
  DMDAVecRestoreArrayRead(da_T,local_Tsource,&Tsource);
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
  PetscScalar    **Tarray,**Tsource,**f,**udot,**rhs_thermal;
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
  f           = FormLocal_thermal(    &info_T , Tarray , f , udot , rhs_thermal , user );

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

  AppCtx         *user = (AppCtx*)ctx;
  DM             pack  = user->pack;
  DM             da_c , da_T;
  DMDALocalInfo  info_c , info_T;
  PetscScalar    **carray,**Tarray,**Tsource,**f_c,**f_T, **udot_c, **udot_T, **eps_2_array, **sigma_array, **rhs_thermal, **rhs_c;
  Vec            local_c, local_T, local_Tsource, local_eps_2, local_sigma, local_Trhs, local_crhs, U_c , U_T , Udot_c , Udot_T , F_c , F_T;
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
  DMGetLocalVector( da_T , &local_crhs);
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
  DMDAVecGetArray(     da_T , local_crhs , &rhs_c );
  DMDAVecGetArrayRead( da_T , local_udotT , &udot_T );
  DMDAVecGetArray(     da_T , local_T , &Tarray );
  DMDAVecGetArrayRead( da_T , local_Tsource , &Tsource );
  DMDAVecGetArray(     da_T , local_Trhs , &rhs_thermal );
  DMDAVecGetArray(     da_T , local_fT , &f_T );
  
  /* Compute function over the locally owned part of the grid */
  compute_eps2_and_sigma_from_temperature( user , U );
  rhs_c       = FormLocalRHS_CH(      &info_c , carray , rhs_c , eps_2_array , sigma_array , user );
  f_c         = FormLocal_CH(         &info_c , carray , f_c , udot_c , rhs_c , user );
  rhs_thermal = FormLocalRHS_thermal( &info_T , Tarray , rhs_thermal , Tsource , user );
  f_T         = FormLocal_thermal(    &info_T , Tarray , f_T , udot_T , rhs_thermal , user );
  
  /* Restore vectors */
  DMDAVecRestoreArray(da_c,local_c,&carray);
  DMDAVecRestoreArrayRead(da_c,local_eps_2,&eps_2_array);
  DMDAVecRestoreArrayRead(da_c,local_sigma,&sigma_array);
  DMDAVecRestoreArrayRead(da_c,local_udotC,&udot_c);
  DMDAVecRestoreArray(da_T,local_T,&Tarray);
  DMDAVecRestoreArrayRead(da_T,local_Tsource,&Tsource);
  DMDAVecRestoreArray(da_T,local_Trhs,&rhs_thermal);
  DMDAVecRestoreArray(da_c,local_crhs,&rhs_c);
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
  DMRestoreLocalVector(da_c,&local_crhs);
  DMRestoreLocalVector(da_c,&local_eps_2);
  DMRestoreLocalVector(da_c,&local_sigma);

  DMRestoreLocalVector(da_c,&local_udotC);
  DMRestoreLocalVector(da_c,&local_fC);
  DMRestoreLocalVector(da_T,&local_udotT);
  DMRestoreLocalVector(da_T,&local_fT);  
  
  PetscFunctionReturn(0);
  
}

// FUNCTIONS FOR THE SPLIT-CH SOLVER

PetscErrorCode FormRHS_CH_split(TS ts,PetscReal t,Vec U,Vec F,void *ctx) {

  // Computes F = RHSfunction

  AppCtx         *user = (AppCtx*)ctx;
  DM              pack = (DM)user->pack;
  DM              da_c , da_T;
  DMDALocalInfo   info_c;
  PetscScalar     **carray,**eps2,**sigma,**rhs_c;
  Vec             local_c,local_eps2,local_sigma,local_rhs_c;
  DMDALocalInfo   info_T;
  PetscScalar     **Tarray,**Tsource,**rhs_phi;
  Vec             local_T,local_Tsource,local_rhs_T;
  Vec             U_c , U_T , F_c , F_T;
  
  PetscFunctionBeginUser;

  // Get composite stuff
  DMCompositeGetEntries( pack , &da_c , &da_T );
  DMCompositeGetAccess(  pack , U     , &U_c , &U_T );
  DMCompositeGetAccess(  pack , F     , &F_c , &F_T );

  // Get CH data types
  DMGetLocalVector( da_c , &local_c );
  DMGetLocalVector( da_c , &local_eps2 );
  DMGetLocalVector( da_c , &local_rhs_c );
  
  DMDAGetLocalInfo( da_c , &info_c );
  
  DMGlobalToLocalBegin( da_c , U_c , INSERT_VALUES , local_c );
  DMGlobalToLocalEnd(   da_c , U_c , INSERT_VALUES , local_c );
  DMGlobalToLocalBegin( da_c , F_c , INSERT_VALUES , local_rhs_c );
  DMGlobalToLocalEnd(   da_c , F_c , INSERT_VALUES , local_rhs_c );
  DMGlobalToLocalBegin( da_c , user->eps_2 , INSERT_VALUES , local_eps2 );
  DMGlobalToLocalEnd(   da_c , user->eps_2 , INSERT_VALUES , local_eps2 );
  
  DMDAVecGetArrayRead( da_c , local_c , &carray );
  DMDAVecGetArrayRead( da_c , local_eps2 , &eps2 );
  DMDAVecGetArray(     da_c , local_rhs_c , &rhs_c );

  // Get phi data types
  DMGetLocalVector( da_T , &local_T );
  DMGetLocalVector( da_T , &local_rhs_T );
  DMGetLocalVector( da_T , &local_sigma );
  
  DMDAGetLocalInfo( da_T , &info_T );
  
  DMGlobalToLocalBegin( da_T , U_T , INSERT_VALUES , local_T );
  DMGlobalToLocalEnd(   da_T , U_T , INSERT_VALUES , local_T );
  DMGlobalToLocalBegin( da_T , F_T , INSERT_VALUES , local_rhs_T );
  DMGlobalToLocalEnd(   da_T , F_T , INSERT_VALUES , local_rhs_T );
  DMGlobalToLocalBegin( da_T , user->sigma , INSERT_VALUES , local_sigma );
  DMGlobalToLocalEnd(   da_T , user->sigma , INSERT_VALUES , local_sigma );

  DMDAVecGetArrayRead( da_T , local_T     , &Tarray );
  DMDAVecGetArrayRead( da_T , local_sigma , &sigma );
  DMDAVecGetArray(     da_T , local_rhs_T , &rhs_phi );
  
  /* Compute function over the locally owned part of the grid */
  rhs_c       = FormLocalRHS_CH_split_c(   &info_c , carray , Tarray, rhs_c , sigma , user );
  rhs_phi     = FormLocalRHS_CH_split_phi( &info_T , carray , rhs_phi , eps2 , user );
  
  /* Restore vectors */
  DMDAVecRestoreArrayRead( da_c , local_c , &carray );
  DMDAVecRestoreArrayRead( da_c , local_eps2 , &eps2 );
  DMDAVecRestoreArray(     da_c , local_rhs_c , &rhs_c );
  DMDAVecRestoreArrayRead( da_T , local_sigma , &sigma );
  DMDAVecRestoreArrayRead( da_T , local_T , &Tarray );
  DMDAVecRestoreArray(     da_T , local_rhs_T , &rhs_phi );
  
  DMLocalToGlobalBegin( da_c , local_rhs_c , INSERT_VALUES , F_c );
  DMLocalToGlobalEnd(   da_c , local_rhs_c , INSERT_VALUES , F_c );
  DMLocalToGlobalBegin( da_T , local_rhs_T , INSERT_VALUES , F_T );
  DMLocalToGlobalEnd(   da_T , local_rhs_T , INSERT_VALUES , F_T );

  DMCompositeRestoreAccess( pack , F , &F_c , &F_T );
  DMCompositeRestoreAccess( pack , U , &U_c , &U_T );

  DMRestoreLocalVector( da_c , &local_c );
  DMRestoreLocalVector( da_c , &local_eps2 );
  DMRestoreLocalVector( da_c , &local_rhs_c );
  DMRestoreLocalVector( da_T , &local_sigma );
  DMRestoreLocalVector( da_T , &local_T );
  DMRestoreLocalVector( da_T , &local_rhs_T );
  
  PetscFunctionReturn(0);
  
}

PetscScalar** FormLocalRHS_CH_split_c( DMDALocalInfo *info ,
                                       PetscScalar **uarray ,
                                       PetscScalar **phiarray ,
                                       PetscScalar **rhs ,
                                       PetscScalar **sigma_array ,
                                       AppCtx *user ) {

  // Function to evaluate RHS of CH dynamics for a given local process
  // Output: rhs array, containing RHS evaluations on local process
  
  // NOTE: these CH eqns are dimensionless with domain length scale = 1. Physical domain size shows up in L_omega.
  PetscScalar hx = 1.0 / (PetscReal)(info->mx-1);
  PetscScalar sx = 1.0/(hx*hx);
  PetscScalar hy = 1.0 / (PetscReal)(info->my-1);
  PetscScalar sy = 1.0/(hy*hy);

  // Set boundary node function depending on user-defined BCs
  void (*set_boundary_ghost_nodes)( AppCtx* , PetscScalar** , PetscInt , PetscInt , PetscInt , PetscInt );
  if (user->boundary.compare("dirichlet") == 0) {
    set_boundary_ghost_nodes = set_boundary_ghost_nodes_dirichlet_singleframe;
  }
  else {
    // NEED TO IMPLEMENT THIS
  }

  /* Compute function over the locally owned part of the grid */
  for (int j=info->ys; j<info->ys+info->ym; j++) {
    for (int i=info->xs; i<info->xs+info->xm; i++) {

      (*set_boundary_ghost_nodes)( user , uarray , info->mx , info->my , i , j );
      
      // dc/dt = laplacian( phi ) - sigma*(c - m)
      PetscScalar sigma = sigma_array[j][i];
      PetscScalar m     = user->m;
      
      // Term: laplacian( phi )
      PetscScalar l_i     = phiarray[j][i];
      PetscScalar l_im1   = phiarray[j][i-1];
      PetscScalar l_ip1   = phiarray[j][i+1];
      PetscScalar l_jm1   = phiarray[j-1][i];
      PetscScalar l_jp1   = phiarray[j+1][i];

      PetscScalar dxx     = sx * ( l_ip1 + l_im1 - 2.0 * l_i );
      PetscScalar dyy     = sy * ( l_jp1 + l_jm1 - 2.0 * l_i );
      
      rhs[j][i]  = dxx + dyy;
      
      // Term: -sigma*(c - m)
      rhs[j][i] += -sigma * ( uarray[j][i] - m );

    }

  }

  return rhs;

}

PetscScalar** FormLocalRHS_CH_split_phi( DMDALocalInfo *info ,
                                         PetscScalar **uarray ,
                                         PetscScalar **rhs ,
                                         PetscScalar **eps2_array ,
                                         AppCtx *user ) {

  // Function to evaluate RHS of CH dynamics for a given local process
  // Output: rhs array, containing RHS evaluations on local process
  
  // NOTE: these CH eqns are dimensionless with domain length scale = 1. Physical domain size shows up in L_omega.
  PetscScalar hx = 1.0 / (PetscReal)(info->mx-1);
  PetscScalar sx = 1.0/(hx*hx);
  PetscScalar hy = 1.0 / (PetscReal)(info->my-1);
  PetscScalar sy = 1.0/(hy*hy);

  // Set boundary node function depending on user-defined BCs
  void (*set_boundary_ghost_nodes)( AppCtx* , PetscScalar** , PetscInt , PetscInt , PetscInt , PetscInt );
  if (user->boundary.compare("dirichlet") == 0) {
    set_boundary_ghost_nodes = set_boundary_ghost_nodes_dirichlet_singleframe;
  }
  else {
    // NEED TO IMPLEMENT THIS
  }

  /* Compute function over the locally owned part of the grid */
  for (int j=info->ys; j<info->ys+info->ym; j++) {
    for (int i=info->xs; i<info->xs+info->xm; i++) {

      (*set_boundary_ghost_nodes)( user , uarray , info->mx , info->my , i , j );
      
      // dphi/dt = -eps2 * laplacian( u ) + ( u^3 - u )
      
      // Term: -eps2 * laplacian( u )
      PetscScalar l_i     = uarray[j][i];
      PetscScalar l_im1   = uarray[j][i-1];
      PetscScalar l_ip1   = uarray[j][i+1];
      PetscScalar l_jm1   = uarray[j-1][i];
      PetscScalar l_jp1   = uarray[j+1][i];

      PetscScalar dxx     = sx * ( l_ip1 + l_im1 - 2.0 * l_i );
      PetscScalar dyy     = sy * ( l_jp1 + l_jm1 - 2.0 * l_i );
      
      rhs[j][i]  = -eps2_array[j][i] * ( dxx + dyy );
      
      // Term: ( u^3 - u )
      rhs[j][i] += uarray[j][i] * uarray[j][i] * uarray[j][i] - uarray[j][i];
      
    }

  }

  return rhs;

}

PetscErrorCode FormIFunction_CH_split(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx) {

  // Computes residual F = Udot - RHSFunction

  AppCtx         *user = (AppCtx*)ctx;
  DM             pack  = user->pack;
  DM             da_c , da_phi , da_T;
  DMDALocalInfo  info_c , info_phi , info_T;
  PetscScalar    **carray,**phiarray,**Tarray,**Tsource,**f_c,**f_phi,**f_T, **udot_c, **udot_phi, **udot_T, **eps_2_array, **sigma_array, **rhs_thermal, **rhs_c, **rhs_phi;
  Vec            local_c, local_phi , local_T, local_Tsource, local_eps_2, local_sigma, local_Trhs, local_crhs, local_phirhs , U_c , U_phi , U_T , Udot_c , Udot_phi , Udot_T , F_c , F_phi , F_T;
  Vec            local_udotC , local_udotPhi , local_udotT , local_fC , local_fPhi , local_fT;
  
  PetscFunctionBeginUser;

  DMCompositeGetEntries( pack , &da_c    , &da_phi   , &da_T );
  DMCompositeGetAccess(  pack , U        , &U_c      , &U_phi    , &U_T );
  DMCompositeGetAccess(  pack , Udot     , &Udot_c   , &Udot_phi , &Udot_T );
  DMCompositeGetAccess(  pack , F        , &F_c      , &F_phi    , &F_T );
  
  DMGetLocalVector( da_c   , &local_c);
  DMGetLocalVector( da_phi , &local_phi);
  DMGetLocalVector( da_T   , &local_T);
  DMGetLocalVector( da_c   , &local_crhs);
  DMGetLocalVector( da_phi , &local_phirhs);
  DMGetLocalVector( da_T   , &local_Trhs);
  DMGetLocalVector( da_phi , &local_eps_2);
  DMGetLocalVector( da_c   , &local_sigma);
  DMGetLocalVector( da_c   , &local_udotC);
  DMGetLocalVector( da_c   , &local_fC);
  DMGetLocalVector( da_phi , &local_udotPhi);
  DMGetLocalVector( da_phi , &local_fPhi);
  DMGetLocalVector( da_T   , &local_udotT);
  DMGetLocalVector( da_T   , &local_fT);  
  DMGetLocalVector( da_T   , &local_Tsource );

  DMDAGetLocalInfo( da_c   , &info_c );
  DMDAGetLocalInfo( da_phi , &info_phi );
  DMDAGetLocalInfo( da_T   , &info_T );
  
  DMGlobalToLocalBegin( da_c , U_c , INSERT_VALUES , local_c );
  DMGlobalToLocalEnd(   da_c , U_c , INSERT_VALUES , local_c );
  DMGlobalToLocalBegin( da_c , Udot_c , INSERT_VALUES , local_udotC );
  DMGlobalToLocalEnd(   da_c , Udot_c , INSERT_VALUES , local_udotC );
  DMGlobalToLocalBegin( da_c , F_c , INSERT_VALUES , local_fC );
  DMGlobalToLocalEnd(   da_c , F_c , INSERT_VALUES , local_fC );
  DMGlobalToLocalBegin( da_c , user->sigma , INSERT_VALUES , local_sigma );
  DMGlobalToLocalEnd(   da_c , user->sigma , INSERT_VALUES , local_sigma );
  DMGlobalToLocalBegin( da_T , U_T , INSERT_VALUES , local_T );
  DMGlobalToLocalEnd(   da_T , U_T , INSERT_VALUES , local_T );
  DMGlobalToLocalBegin( da_T , Udot_T , INSERT_VALUES , local_udotT );
  DMGlobalToLocalEnd(   da_T , Udot_T , INSERT_VALUES , local_udotT );
  DMGlobalToLocalBegin( da_T , F_T , INSERT_VALUES , local_fT );
  DMGlobalToLocalEnd(   da_T , F_T , INSERT_VALUES , local_fT );
  DMGlobalToLocalBegin( da_T , user->temperature_source , INSERT_VALUES , local_Tsource );
  DMGlobalToLocalEnd(   da_T , user->temperature_source , INSERT_VALUES , local_Tsource );  
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
  DMDAVecGetArray(     da_phi , local_phi   , &phiarray );
  DMDAVecGetArrayRead( da_phi , local_eps_2 , &eps_2_array );
  DMDAVecGetArrayRead( da_phi , local_udotPhi , &udot_phi );
  DMDAVecGetArray(     da_phi , local_fPhi , &f_phi );
  DMDAVecGetArray(     da_phi , local_phirhs , &rhs_phi );
  DMDAVecGetArray(     da_T , local_crhs , &rhs_c );
  DMDAVecGetArrayRead( da_T , local_udotT , &udot_T );
  DMDAVecGetArray(     da_T , local_T , &Tarray );
  DMDAVecGetArray(     da_T , local_Trhs , &rhs_thermal );
  DMDAVecGetArray(     da_T , local_fT , &f_T );
  DMDAVecGetArrayRead( da_T , local_Tsource , &Tsource );
  
  /* Compute function over the locally owned part of the grid */
  //compute_eps2_and_sigma_from_temperature( user , U );
  rhs_c       = FormLocalRHS_CH_split_c(    &info_c   , carray   , phiarray , rhs_c , sigma_array , user );
  f_c         = FormLocal_CH(               &info_c   , carray   , f_c      , udot_c , rhs_c , user );
  rhs_phi     = FormLocalRHS_CH_split_phi(  &info_phi , carray   , rhs_phi  , eps_2_array , user );
  f_phi       = FormLocal_CH(               &info_phi , phiarray , f_phi    , phiarray    , rhs_phi , user );
  rhs_thermal = FormLocalRHS_thermal(       &info_T   , Tarray   , rhs_thermal , Tsource , user );
  f_T         = FormLocal_CH(               &info_T   , Tarray   , f_T , udot_T , rhs_thermal , user );
  
  /* Restore vectors */
  DMDAVecRestoreArray(da_c,local_c,&carray);
  DMDAVecRestoreArrayRead(da_c,local_sigma,&sigma_array);
  DMDAVecRestoreArrayRead(da_c,local_udotC,&udot_c);
  DMDAVecRestoreArray(da_c,local_crhs,&rhs_c);
  DMDAVecRestoreArrayRead(da_phi,local_eps_2,&eps_2_array);
  DMDAVecRestoreArray(da_phi,local_phi,&phiarray);
  DMDAVecRestoreArrayRead(da_phi,local_udotPhi,&udot_phi);
  DMDAVecRestoreArray(da_phi,local_phirhs,&rhs_phi);
  DMDAVecRestoreArray(da_T,local_T,&Tarray);
  DMDAVecRestoreArray(da_T,local_Trhs,&rhs_thermal);
  DMDAVecRestoreArrayRead(da_T,local_udotT,&udot_T);

  DMDAVecRestoreArray(da_c  ,local_fC  ,&f_c);
  DMDAVecRestoreArray(da_phi,local_fPhi,&f_phi);
  DMDAVecRestoreArray(da_T  ,local_fT  ,&f_T);

  DMLocalToGlobalBegin( da_c   , local_fC   , INSERT_VALUES , F_c );
  DMLocalToGlobalEnd(   da_c   , local_fC   , INSERT_VALUES , F_c );
  DMLocalToGlobalBegin( da_phi , local_fPhi , INSERT_VALUES , F_phi );
  DMLocalToGlobalEnd(   da_phi , local_fPhi , INSERT_VALUES , F_phi );
  DMLocalToGlobalBegin( da_T   , local_fT   , INSERT_VALUES , F_T );
  DMLocalToGlobalEnd(   da_T   , local_fT   , INSERT_VALUES , F_T );

  DMCompositeRestoreAccess(  pack , U        , &U_c    , &U_phi    , &U_T );
  DMCompositeRestoreAccess(  pack , Udot     , &Udot_c , &Udot_phi , &Udot_T );
  DMCompositeRestoreAccess(  pack , F        , &F_c    , &F_phi    , &F_T );
  
  DMRestoreLocalVector(da_c,&local_c);
  DMRestoreLocalVector(da_c,&local_crhs);
  DMRestoreLocalVector(da_c,&local_sigma);
  DMRestoreLocalVector(da_phi,&local_phi);
  DMRestoreLocalVector(da_phi,&local_phirhs);
  DMRestoreLocalVector(da_phi,&local_eps_2);
  DMRestoreLocalVector(da_T,&local_T);
  DMRestoreLocalVector(da_T,&local_Trhs);

  DMRestoreLocalVector(da_c,&local_udotC);
  DMRestoreLocalVector(da_c,&local_fC);
  DMRestoreLocalVector(da_phi,&local_udotPhi);
  DMRestoreLocalVector(da_phi,&local_fPhi);
  DMRestoreLocalVector(da_T,&local_udotT);
  DMRestoreLocalVector(da_T,&local_fT);  
  
  PetscFunctionReturn(0);
  
}
