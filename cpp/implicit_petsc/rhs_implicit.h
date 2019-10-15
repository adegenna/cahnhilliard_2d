#ifndef __RHS_IMPLICIT_H__
#define __RHS_IMPLICIT_H__

#include "utils_ch_implicit.h"

PetscErrorCode FormLocal_CH( DMDALocalInfo *info ,
                             PetscScalar **uarray ,
                             PetscScalar **eps_2_array ,
                             PetscScalar **sigma_array ,
                             PetscScalar **f , 
                             PetscScalar** udot ,
                             AppCtx *ctx );

PetscScalar** FormLocalRHS_thermal( DMDALocalInfo *info ,
                                    PetscScalar **Tarray ,
                                    PetscScalar **rhs ,
                                    PetscScalar **Tsource ,
                                    AppCtx *user );

PetscScalar** FormLocal_thermal( DMDALocalInfo *info ,
				 PetscScalar **Tarray ,
				 PetscScalar **f , 
				 PetscScalar** udot ,
				 AppCtx *ctx );

PetscErrorCode FormRHS_thermal(TS ts,PetscReal t,Vec U,Vec F,void *ctx);

PetscErrorCode FormIFunction_CH( TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx );

PetscErrorCode FormIFunction_CH_coupled(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx);

PetscErrorCode FormIFunction_thermal(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx);

#endif
