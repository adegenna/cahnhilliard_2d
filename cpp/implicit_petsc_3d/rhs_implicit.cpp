#include "rhs_implicit.h"
#include "boundary_conditions.h"
#include <petscdmcomposite.h>

PetscScalar*** FormLocalResidual_ch( DMDALocalInfo *info ,
                                     PetscScalar ***uarray ,
                                     PetscScalar ***u_optional ,
                                     PetscScalar ***f , 
                                     PetscScalar ***udot ,
                                     PetscScalar ***rhs ,
                                     AppCtx *user ) {
  
  // Apply correction to boundary fluxes and set residuals
  for (int k=info->zs; k<info->zs+info->zm; k++) {
    for (int j=info->ys; j<info->ys+info->ym; j++) {
      for (int i=info->xs; i<info->xs+info->xm; i++) {

        f[k][j][i] = user->residualFunction_ch( uarray , u_optional , rhs[k][j][i] , udot[k][j][i] , info->mx , info->my , info->mz , i , j , k );
      
      }
    }
  }
  
  return f;

}

PetscScalar*** FormLocalResidual_thermal( DMDALocalInfo *info ,
                                          PetscScalar ***uarray ,
                                          PetscScalar ***u_optional ,
                                          PetscScalar ***f , 
                                          PetscScalar ***udot ,
                                          PetscScalar ***rhs ,
                                          AppCtx *user ) {
  
  // Apply correction to boundary fluxes and set residuals
  for (int k=info->zs; k<info->zs+info->zm; k++) {
    for (int j=info->ys; j<info->ys+info->ym; j++) {
      for (int i=info->xs; i<info->xs+info->xm; i++) {

        f[k][j][i] = user->residualFunction_thermal( uarray , u_optional , rhs[k][j][i] , udot[k][j][i] , info->mx , info->my , info->mz , i , j , k );
      
      }
    }
  }
  
  return f;

}

PetscScalar*** FormLocalRHS_thermal( DMDALocalInfo *info ,
                                     PetscScalar ***Tarray ,
                                     PetscScalar ***rhs ,
                                     PetscScalar ***Tsource ,
                                     AppCtx *user ) {
  
  // Function to evaluate RHS of thermal dynamics for a given local process
  // Output: rhs array, containing RHS evaluations on local process

  // NOTE: these CH eqns are dimensionless with domain length scale = 1. Physical domain size shows up in L_omega.
  PetscScalar hx = 1.0 / (PetscReal)(info->mx-1);
  PetscScalar sx = 1.0 / (hx*hx);
  PetscScalar hy = 1.0 / (PetscReal)(info->my-1);
  PetscScalar sy = 1.0 / (hy*hy);
  PetscScalar hz = 1.0 / (PetscReal)(info->mz-1);
  PetscScalar sz = 1.0 / (hz*hz);
  
  /* Compute function over the locally owned part of the grid */
  for (int k=info->zs; k<info->zs+info->zm; k++) {
    for (int j=info->ys; j<info->ys+info->ym; j++) {
      for (int i=info->xs; i<info->xs+info->xm; i++) {
        
        // dT/dt = D_T * laplacian( T ) + S
      
        // Term: D_T * laplacian( T )
        PetscScalar dxx     = sx * ( Tarray[k][j][i+1] + Tarray[k][j][i-1] - 2.0 * Tarray[k][j][i] );
        PetscScalar dyy     = sy * ( Tarray[k][j+1][i] + Tarray[k][j-1][i] - 2.0 * Tarray[k][j][i] );
        PetscScalar dzz     = sz * ( Tarray[k+1][j][i] + Tarray[k-1][j][i] - 2.0 * Tarray[k][j][i] );
      
        rhs[k][j][i]        = user->D_T * ( dxx + dyy + dzz ) + Tsource[k][j][i];
      
      }
    }
  }

  return rhs;
  
}

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

  /* Compute function over the locally owned part of the grid */
  for (int k=info->zs; k<info->zs+info->zm; k++) {
    for (int j=info->ys; j<info->ys+info->ym; j++) {
      for (int i=info->xs; i<info->xs+info->xm; i++) {
        
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

  /* Compute function over the locally owned part of the grid */
  for (int k=info->zs; k<info->zs+info->zm; k++) {
    for (int j=info->ys; j<info->ys+info->ym; j++) {
      for (int i=info->xs; i<info->xs+info->xm; i++) {
        
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
  PetscScalar    ***carray,***phiarray, ***f_c,***f_phi, ***udot_c,***udot_phi, ***eps_2_array,***sigma_array, ***rhs_c, ***rhs_phi, ***ch_bc_array;
  Vec            local_c,local_phi, local_eps_2,local_sigma, local_crhs,local_phirhs, U_c,U_phi, Udot_c,Udot_phi, F_c , F_phi;
  Vec            local_udotC,local_udotPhi, local_fC,local_fPhi,local_ch_bc_array;
  
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
  DMGetLocalVector( da_c   , &local_ch_bc_array);

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
  DMGlobalToLocalBegin( da_c , user->dirichlet_bc_ch_array , INSERT_VALUES , local_ch_bc_array );
  DMGlobalToLocalEnd(   da_c , user->dirichlet_bc_ch_array , INSERT_VALUES , local_ch_bc_array );
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
  DMDAVecGetArrayRead( da_c , local_ch_bc_array , &ch_bc_array );
  DMDAVecGetArray(     da_phi , local_phi   , &phiarray );
  DMDAVecGetArrayRead( da_phi , local_eps_2 , &eps_2_array );
  DMDAVecGetArrayRead( da_phi , local_udotPhi , &udot_phi );
  DMDAVecGetArray(     da_phi , local_fPhi , &f_phi );
  DMDAVecGetArray(     da_phi , local_phirhs , &rhs_phi );
  
  /* Compute function over the locally owned part of the grid */
  rhs_c       = FormLocalRHS_CH_split_c(    &info_c   , carray   , phiarray    , rhs_c , sigma_array , user );
  f_c         = FormLocalResidual_ch(       &info_c   , carray   , ch_bc_array , f_c      , udot_c , rhs_c , user );
  rhs_phi     = FormLocalRHS_CH_split_phi(  &info_phi , carray   , rhs_phi     , eps_2_array , user );
  f_phi       = FormLocalResidual_ch(       &info_phi , phiarray , NULL        , f_phi    , phiarray    , rhs_phi , user );
  
  /* Restore vectors */
  DMDAVecRestoreArray(da_c,local_c,&carray);
  DMDAVecRestoreArrayRead(da_c,local_sigma,&sigma_array);
  DMDAVecRestoreArrayRead(da_c,local_udotC,&udot_c);
  DMDAVecRestoreArray(da_c,local_crhs,&rhs_c);
  DMDAVecRestoreArrayRead(da_c,local_ch_bc_array,&ch_bc_array);
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
  DMRestoreLocalVector(da_c,&local_ch_bc_array);
  DMRestoreLocalVector(da_phi,&local_phi);
  DMRestoreLocalVector(da_phi,&local_phirhs);
  DMRestoreLocalVector(da_phi,&local_eps_2);

  DMRestoreLocalVector(da_c,&local_udotC);
  DMRestoreLocalVector(da_c,&local_fC);
  DMRestoreLocalVector(da_phi,&local_udotPhi);
  DMRestoreLocalVector(da_phi,&local_fPhi);
  
  PetscFunctionReturn(0);
  
}



PetscErrorCode FormIFunction_CH_split_thermal(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx) {

  // Computes residual F = Udot - RHSFunction

  AppCtx         *user = (AppCtx*)ctx;
  DM             pack  = user->pack;
  DM             da_c , da_phi , da_T;
  DMDALocalInfo  info_c , info_phi , info_T;
  PetscScalar    ***carray,***phiarray,***Tarray, ***f_c,***f_phi,***f_T, ***udot_c,***udot_phi,***udot_T, ***eps_2_array,***sigma_array, ***rhs_c,***rhs_phi,***rhs_T, ***Tsource , ***thermal_bc_array,***ch_bc_array;
  Vec            local_c,local_phi,local_T, local_eps_2,local_sigma, local_crhs,local_phirhs,local_Trhs, U_c,U_phi,U_T, Udot_c,Udot_phi,Udot_T, F_c,F_phi,F_T;
  Vec            local_udotC,local_udotPhi,local_udotT, local_fC,local_fPhi,local_fT,local_Tsource,local_thermal_bc_array,local_ch_bc_array;
  
  PetscFunctionBeginUser;

  DMCompositeGetEntries( pack , &da_c    , &da_phi , &da_T );
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
  DMGetLocalVector( da_T   , &local_thermal_bc_array );
  DMGetLocalVector( da_c   , &local_ch_bc_array );

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
  DMGlobalToLocalBegin( da_phi , U_phi , INSERT_VALUES , local_phi );
  DMGlobalToLocalEnd(   da_phi , U_phi , INSERT_VALUES , local_phi );
  DMGlobalToLocalBegin( da_phi , Udot_phi , INSERT_VALUES , local_udotPhi );
  DMGlobalToLocalEnd(   da_phi , Udot_phi , INSERT_VALUES , local_udotPhi );
  DMGlobalToLocalBegin( da_phi , F_phi , INSERT_VALUES , local_fPhi );
  DMGlobalToLocalEnd(   da_phi , F_phi , INSERT_VALUES , local_fPhi );
  DMGlobalToLocalBegin( da_phi , user->eps_2 , INSERT_VALUES , local_eps_2 );
  DMGlobalToLocalEnd(   da_phi , user->eps_2 , INSERT_VALUES , local_eps_2 );
  DMGlobalToLocalBegin( da_T , U_T , INSERT_VALUES , local_T );
  DMGlobalToLocalEnd(   da_T , U_T , INSERT_VALUES , local_T );
  DMGlobalToLocalBegin( da_T , Udot_T , INSERT_VALUES , local_udotT );
  DMGlobalToLocalEnd(   da_T , Udot_T , INSERT_VALUES , local_udotT );
  DMGlobalToLocalBegin( da_T , F_T , INSERT_VALUES , local_fT );
  DMGlobalToLocalEnd(   da_T , F_T , INSERT_VALUES , local_fT );
  DMGlobalToLocalBegin( da_T , user->temperature_source , INSERT_VALUES , local_Tsource );
  DMGlobalToLocalEnd(   da_T , user->temperature_source , INSERT_VALUES , local_Tsource );
  DMGlobalToLocalBegin( da_T , user->dirichlet_bc_thermal_array , INSERT_VALUES , local_thermal_bc_array );
  DMGlobalToLocalEnd(   da_T , user->dirichlet_bc_thermal_array , INSERT_VALUES , local_thermal_bc_array );
  DMGlobalToLocalBegin( da_c , user->dirichlet_bc_ch_array , INSERT_VALUES , local_ch_bc_array );
  DMGlobalToLocalEnd(   da_c , user->dirichlet_bc_ch_array , INSERT_VALUES , local_ch_bc_array );
  
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
  DMDAVecGetArray(     da_T , local_Trhs    , &rhs_T );
  DMDAVecGetArrayRead( da_T , local_udotT   , &udot_T );
  DMDAVecGetArray(     da_T , local_T       , &Tarray );
  DMDAVecGetArray(     da_T , local_Trhs    , &rhs_T );
  DMDAVecGetArray(     da_T , local_fT      , &f_T );
  DMDAVecGetArrayRead( da_T , local_Tsource , &Tsource );
  DMDAVecGetArrayRead( da_T , local_thermal_bc_array , &thermal_bc_array );
  DMDAVecGetArrayRead( da_c , local_ch_bc_array , &ch_bc_array );

  /* Compute function over the locally owned part of the grid */
  rhs_c       = FormLocalRHS_CH_split_c(    &info_c   , carray   , phiarray         , rhs_c       , sigma_array , user );
  f_c         = FormLocalResidual_ch(       &info_c   , carray   , ch_bc_array      , f_c          , udot_c      , rhs_c       , user );
  rhs_phi     = FormLocalRHS_CH_split_phi(  &info_phi , carray   , rhs_phi          , eps_2_array , user );
  f_phi       = FormLocalResidual_ch(       &info_phi , phiarray , NULL             , f_phi        , phiarray    , rhs_phi     , user );
  rhs_T       = FormLocalRHS_thermal(       &info_T   , Tarray   , rhs_T            , Tsource     , user );
  f_T         = FormLocalResidual_thermal(  &info_T   , Tarray   , thermal_bc_array , f_T          , udot_T      , rhs_T       , user );

  /* Restore vectors */
  DMDAVecRestoreArray(     da_c   , local_c       , &carray);
  DMDAVecRestoreArrayRead( da_c   , local_sigma   , &sigma_array);
  DMDAVecRestoreArrayRead( da_c   , local_udotC   , &udot_c);
  DMDAVecRestoreArray(     da_c   , local_crhs    , &rhs_c);
  DMDAVecRestoreArrayRead( da_phi , local_eps_2   , &eps_2_array);
  DMDAVecRestoreArray(     da_phi , local_phi     , &phiarray);
  DMDAVecRestoreArrayRead( da_phi , local_udotPhi , &udot_phi);
  DMDAVecRestoreArray(     da_phi , local_phirhs  , &rhs_phi);
  DMDAVecRestoreArray(     da_T   , local_T       , &Tarray);
  DMDAVecRestoreArray(     da_T   , local_Trhs    , &rhs_T);
  DMDAVecRestoreArrayRead( da_T   , local_udotT   , &udot_T);
  DMDAVecRestoreArrayRead( da_T   , local_Tsource , &Tsource);
  DMDAVecRestoreArrayRead( da_T   , local_thermal_bc_array , &thermal_bc_array );
  DMDAVecRestoreArrayRead( da_T   , local_ch_bc_array , &ch_bc_array );

  DMDAVecRestoreArray(da_c   , local_fC   , &f_c);
  DMDAVecRestoreArray(da_phi , local_fPhi , &f_phi);
  DMDAVecRestoreArray(da_T   , local_fT   , &f_T);

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
  DMRestoreLocalVector(da_c,&local_udotC);
  DMRestoreLocalVector(da_c,&local_fC);
  DMRestoreLocalVector(da_phi,&local_phi);
  DMRestoreLocalVector(da_phi,&local_phirhs);
  DMRestoreLocalVector(da_phi,&local_eps_2);
  DMRestoreLocalVector(da_phi,&local_udotPhi);
  DMRestoreLocalVector(da_phi,&local_fPhi);
  DMRestoreLocalVector(da_T,&local_T);
  DMRestoreLocalVector(da_T,&local_Trhs);
  DMRestoreLocalVector(da_T,&local_udotT);
  DMRestoreLocalVector(da_T,&local_fT);  
  DMRestoreLocalVector(da_T,&local_Tsource );
  DMRestoreLocalVector(da_T,&local_thermal_bc_array );
  DMRestoreLocalVector(da_c,&local_ch_bc_array );

  PetscFunctionReturn(0);
  
}
