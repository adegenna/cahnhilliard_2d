#ifndef __BOUNDARY_CONDITIONS_H__
#define __BOUNDARY_CONDITIONS_H__

#include "utils_ch_implicit.h"

ThirteenPointStencil apply_dirichlet_bc( AppCtx* user ,
					 PetscReal** uarray ,
					 PetscInt Mx , PetscInt My ,
					 PetscInt i , PetscInt j );



#endif
