#ifndef __RHS_EXPLICIT_H__
#define __RHS_EXPLICIT_H__

#include "utils_ch_explicit.h"

PetscScalar** FormLocalRHS_CH( DMDALocalInfo *info ,
                               PetscScalar **uarray ,
                               PetscScalar **rhs ,
                               PetscScalar **eps_2_array ,
                               PetscScalar **sigma_array ,
                               AppCtx *user );

PetscScalar** FormLocalRHS_thermal( DMDALocalInfo *info ,
                                    PetscScalar **Tarray ,
                                    PetscScalar **rhs ,
                                    PetscScalar **Tsource ,
                                    AppCtx *user );

PetscScalar** set_boundary_values( DMDALocalInfo *info ,
                                   PetscScalar **uarray ,
                                   PetscScalar **u_optional ,
				   AppCtx *user );

PetscErrorCode FormRHS_CH_coupled(TS ts,PetscReal t,Vec U,Vec F,void *ctx);

PetscErrorCode FormRHS_CH(TS ts,PetscReal t,Vec U,Vec F,void *ctx);

PetscErrorCode FormRHS_thermal(TS ts,PetscReal t,Vec U,Vec F,void *ctx);

#endif
