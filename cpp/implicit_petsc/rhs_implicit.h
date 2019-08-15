#ifndef __RHS_IMPLICIT_H__
#define __RHS_IMPLICIT_H__

#include "utils_ch_implicit.h"

PetscErrorCode compute_rhs( TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx );

PetscErrorCode FormIFunction( TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx );


#endif
