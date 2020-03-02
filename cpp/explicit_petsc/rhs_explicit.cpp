#include <petscdmcomposite.h>
#include "rhs_explicit.h"
#include "boundary_conditions.h"
#include "temperature_dependence.h"

PetscScalar** set_boundary_values( DMDALocalInfo *info ,
                                   PetscScalar **uarray ,
                                   PetscScalar **u_optional ,
                                   AppCtx *user ) {
  
  // Apply correction to boundary values
  for (int j=info->ys; j<info->ys+info->ym; j++) {
    for (int i=info->xs; i<info->xs+info->xm; i++) {

      uarray[j][i] = user->residualFunction_ch( uarray , info->mx , info->my , i , j );

    }
  }

  return uarray;
  
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

  /* Compute function over the locally owned part of the grid */
  for (int j=info->ys; j<info->ys+info->ym; j++) {
    for (int i=info->xs; i<info->xs+info->xm; i++) {
      
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

      // Boundary conditions
      if ( i <= 1 || j <= 1 || i >= (info->mx-2) || j >= (info->my-2) )
        rhs[j][i] = 0.0;

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
  PetscScalar     **carray,**eps2,**sigma,**rhs_c , **ch_bc_array;
  Vec             local_c,local_eps2,local_sigma,local_rhs , local_ch_bc_array;
  
  PetscFunctionBeginUser;

  DMCompositeGetEntries( pack , &da_c , &da_T );
  
  DMGetLocalVector( da_c , &local_c );
  DMGetLocalVector( da_c , &local_eps2 );
  DMGetLocalVector( da_c , &local_sigma );
  DMGetLocalVector( da_c , &local_rhs );
  DMGetLocalVector( da_c , &local_ch_bc_array);

  DMDAGetLocalInfo( da_c , &info_c );
  
  DMGlobalToLocalBegin( da_c , U , INSERT_VALUES , local_c );
  DMGlobalToLocalEnd(   da_c , U , INSERT_VALUES , local_c );
  DMGlobalToLocalBegin( da_c , F , INSERT_VALUES , local_rhs );
  DMGlobalToLocalEnd(   da_c , F , INSERT_VALUES , local_rhs );
  DMGlobalToLocalBegin( da_c , user->eps_2 , INSERT_VALUES , local_eps2 );
  DMGlobalToLocalEnd(   da_c , user->eps_2 , INSERT_VALUES , local_eps2 );
  DMGlobalToLocalBegin( da_c , user->sigma , INSERT_VALUES , local_sigma );
  DMGlobalToLocalEnd(   da_c , user->sigma , INSERT_VALUES , local_sigma );
  DMGlobalToLocalBegin( da_c , user->dirichlet_bc_ch_array , INSERT_VALUES , local_ch_bc_array );
  DMGlobalToLocalEnd(   da_c , user->dirichlet_bc_ch_array , INSERT_VALUES , local_ch_bc_array );

  DMDAVecGetArray( da_c , local_c , &carray );
  DMDAVecGetArrayRead( da_c , local_eps2 , &eps2 );
  DMDAVecGetArrayRead( da_c , local_sigma , &sigma );
  DMDAVecGetArray(     da_c , local_rhs , &rhs_c );
  DMDAVecGetArrayRead( da_c , local_ch_bc_array , &ch_bc_array );

  /* Compute function over the locally owned part of the grid */
  rhs_c  = FormLocalRHS_CH( &info_c , carray , rhs_c , eps2 , sigma , user );
  rhs_c = set_boundary_values( &info_c , rhs_c , ch_bc_array , user );
  
  /* Restore vectors */
  DMDAVecRestoreArray( da_c , local_c , &carray );
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
  
  /* Compute function over the locally owned part of the grid */
  for (int j=info->ys; j<info->ys+info->ym; j++) {
    for (int i=info->xs; i<info->xs+info->xm; i++) {
      
      // dT/dt = D_T * laplacian( T ) + S
      
      // Term: D_T * laplacian( T )
      PetscScalar dxx     = sx * ( Tarray[j][i+1] + Tarray[j][i-1] - 2.0 * Tarray[j][i] );
      PetscScalar dyy     = sy * ( Tarray[j+1][i] + Tarray[j-1][i] - 2.0 * Tarray[j][i] );
      
      rhs[j][i]           = user->D_T * ( dxx + dyy ) + Tsource[j][i];

      // Boundary conditions
      if ( i <= 1 || j <= 1 || i >= (info->mx-2) || j >= (info->my-2) )
        rhs[j][i] = 0.0;
      
    }

  }

  return rhs;
  
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
  rhs_c       = set_boundary_values(  &info_c , rhs_c , NULL , user );
  rhs_thermal = FormLocalRHS_thermal( &info_T , Tarray , rhs_thermal , Tsource , user );
  rhs_thermal = set_boundary_values(  &info_T , rhs_thermal , NULL , user );

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
