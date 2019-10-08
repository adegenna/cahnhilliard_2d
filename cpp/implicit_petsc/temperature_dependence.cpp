#include "temperature_dependence.h"
#include "utils_ch_implicit.h"
#include <petscdmcomposite.h>

double convert_temperature_to_flory_huggins( const double T ,
					     const double X_min ,
					     const double X_max ,
					     const double T_min ,
					     const double T_max ) {

  const double dX_dTinv   = ( X_max  - X_min ) / ( 1.0 / T_min - 1.0 / T_max );  
  const double dTinv      = 1.0 / T - 1.0 / T_max;
  double X                = dX_dTinv * dTinv + X_min;

  return X;

}


double compute_eps2_from_chparams( const double X ,
				   const double L_kuhn ,
				   const double m ,
				   const double L_omega ) {

  const double m_scaled   = 0.5 * ( 1.0 - m );
  double eps_2            = L_kuhn * L_kuhn / ( 3.0 * m_scaled * (1.0 - m_scaled) * X * L_omega * L_omega );

  return eps_2;

}

double compute_sigma_from_chparams( const double X ,
				    const double L_kuhn ,
				    const double m ,
				    const double L_omega ,
				    const double N ) {
  
  const double m_scaled   = 0.5 * ( 1.0 - m );
  double sigma            = 36.0 * L_omega * L_omega / ( m_scaled * m_scaled * (1.0 - m_scaled) * (1.0 - m_scaled) * L_kuhn * L_kuhn * X * N * N );  

  return sigma;

}

PetscErrorCode compute_eps2_and_sigma_from_temperature( void *ctx , Vec U ) {

  PetscErrorCode ierr;
  AppCtx         *user = (AppCtx*)ctx;
  DM             pack  = (DM)user->pack;
  DM             da_c , da_T;
  
  PetscInt       i,j,xs,ys,xm,ym,Mx,My;

  PetscScalar    **tarray , **xarray , **eps2array , **sigmaarray;
  Vec            local_temperature , local_X, local_eps2, local_sigma;
  Vec            U_c , U_T;

  PetscFunctionBeginUser;

  DMCompositeGetEntries( pack , &da_c    , &da_T );
  DMCompositeGetAccess(  pack , U        , &U_c , &U_T );
  
  DMGetLocalVector(da_c,&local_X);
  DMGetLocalVector(da_T,&local_temperature);
  DMGetLocalVector(da_c,&local_eps2);
  DMGetLocalVector(da_c,&local_sigma);

  DMDAGetInfo(da_c,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);  

  DMGlobalToLocalBegin( da_T , U_T , INSERT_VALUES , local_temperature );
  DMGlobalToLocalEnd(   da_T , U_T , INSERT_VALUES , local_temperature );
  DMGlobalToLocalBegin( da_c , user->X , INSERT_VALUES , local_X );
  DMGlobalToLocalEnd(   da_c , user->X , INSERT_VALUES , local_X );
  DMGlobalToLocalBegin( da_c , user->eps_2 , INSERT_VALUES , local_eps2 );
  DMGlobalToLocalEnd(   da_c , user->eps_2 , INSERT_VALUES , local_eps2 );
  DMGlobalToLocalBegin( da_c , user->sigma , INSERT_VALUES , local_sigma );
  DMGlobalToLocalEnd(   da_c , user->sigma , INSERT_VALUES , local_sigma );

  /* Get pointers to vector data */
  DMDAVecGetArrayRead( da_T , local_temperature , &tarray );
  DMDAVecGetArray( da_c , local_X , &xarray );
  DMDAVecGetArray( da_c , local_eps2 , &eps2array );
  DMDAVecGetArray( da_c , local_sigma , &sigmaarray );
  
  /* Get local grid boundaries */
  DMDAGetCorners( da_c , &xs , &ys , NULL , &xm , &ym , NULL );
  
  /* Compute function over the locally owned part of the grid */
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      
      xarray[j][i]    = convert_temperature_to_flory_huggins( tarray[j][i] ,
							      user->X_min ,
							      user->X_max ,
							      user->T_min ,
							      user->T_max );
      
      eps2array[j][i] = compute_eps2_from_chparams( xarray[j][i] ,
						    user->L_kuhn ,
						    user->m ,
						    user->L_omega );

      sigmaarray[j][i] = compute_sigma_from_chparams( xarray[j][i] ,
						      user->L_kuhn ,
						      user->m ,
						      user->L_omega ,
						      user->N );

      eps2array[j][i]   = std::min( std::max( eps2array[j][i]  , user->eps2_min )   , user->eps2_max );
      sigmaarray[j][i]  = std::min( std::max( sigmaarray[j][i] , user->sigma_min )  , user->sigma_max );
      
    }
  }

  /* Restore vectors */
  DMDAVecRestoreArrayRead(da_T,local_temperature,&tarray);
  DMDAVecRestoreArray(da_c,local_X,&xarray);
  DMDAVecRestoreArray(da_c,local_eps2,&eps2array);
  DMDAVecRestoreArray(da_c,local_sigma,&sigmaarray);

  DMLocalToGlobalBegin( da_c , local_eps2 , INSERT_VALUES , user->eps_2 );
  DMLocalToGlobalEnd(   da_c , local_eps2 , INSERT_VALUES , user->eps_2 );
  DMLocalToGlobalBegin( da_c , local_sigma , INSERT_VALUES , user->sigma );
  DMLocalToGlobalEnd(   da_c , local_sigma , INSERT_VALUES , user->sigma );
  DMLocalToGlobalBegin( da_c , local_X     , INSERT_VALUES , user->X );
  DMLocalToGlobalEnd(   da_c , local_X     , INSERT_VALUES , user->X );
  
  DMRestoreLocalVector(da_T,&local_temperature);
  DMRestoreLocalVector(da_c,&local_X);
  DMRestoreLocalVector(da_c,&local_eps2);
  DMRestoreLocalVector(da_c,&local_sigma);

  DMCompositeRestoreAccess(  pack , U , &U_c , &U_T );
  
  PetscLogFlops(11.0*ym*xm);
  PetscFunctionReturn(0);
      
}
