#include "temperature_dependence.h"
#include "utils_ch_implicit.h"

double convert_temperature_to_flory_huggins( const double T ,
					     const double X_min ,
					     const double X_max ,
					     const double T_min ,
					     const double T_max ) {

  const double dX_dTinv   = ( X_max  - X_min ) / ( 1.0 / T_min - 1.0 / T_max );  
  const double dTinv      = 1.0 / T - 1.0 / T_max;
  const double X          = dX_dTinv * dTinv + X_min;

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

PetscErrorCode compute_eps2_and_sigma_from_temperature( void *ctx ) {

  PetscErrorCode ierr;
  AppCtx         *user = (AppCtx*)ctx;
  DM             da    = (DM)user->da;

  PetscInt       i,j,xs,ys,xm,ym;

  PetscScalar    **tarray , **xarray , **eps2array , **sigmaarray;
  Vec            local_temperature , local_X;
  
  /*
     Scatter ghost points to local vector,using the 2-step process
        DMGlobalToLocalBegin(),DMGlobalToLocalEnd().
     By placing code between these two statements, computations can be
     done while messages are in transition.
  */
  DMGlobalToLocalBegin( da , user->temperature , INSERT_VALUES , local_temperature );
  DMGlobalToLocalEnd( da , user->temperature , INSERT_VALUES , local_temperature );
  DMGlobalToLocalBegin( da , user->X , INSERT_VALUES , local_X );
  DMGlobalToLocalEnd( da , user->X , INSERT_VALUES , local_X );

  /* Get pointers to vector data */
  DMDAVecGetArrayRead( da , local_temperature , &tarray );
  DMDAVecGetArrayRead( da , local_X , &xarray );

  /* Get local grid boundaries */
  DMDAGetCorners( da , &xs , &ys , NULL , &xm , &ym , NULL );

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
  DMDAVecRestoreArrayRead(da,local_temperature,&tarray);
  DMDAVecRestoreArrayRead(da,local_X,&xarray);
  DMRestoreLocalVector(da,&local_temperature);
  DMRestoreLocalVector(da,&local_X);
  PetscLogFlops(11.0*ym*xm);
  PetscFunctionReturn(0);
      
}
