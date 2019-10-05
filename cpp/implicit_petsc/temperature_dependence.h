#ifndef __TEMPERATURE_DEPENDENCE_H__
#define __TEMPERATURE_DEPENDENCE_H__

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscdmlabel.h>     
#include <petscds.h>
#include <stdio.h>

double convert_temperature_to_flory_huggins( const double T ,
					     const double X_min ,
					     const double X_max ,
					     const double T_min ,
					     const double T_max );

double compute_eps2_from_chparams( const double X ,
				   const double L_kuhn ,
				   const double m ,
				   const double L_omega );

double compute_sigma_from_chparams( const double X ,
				    const double L_kuhn ,
				    const double m ,
				    const double L_omega ,
				    const double N );

PetscErrorCode compute_eps2_and_sigma_from_temperature( void *ctx , Vec U );


#endif
