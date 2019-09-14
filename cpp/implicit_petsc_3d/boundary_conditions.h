#ifndef __BOUNDARY_CONDITIONS_H__
#define __BOUNDARY_CONDITIONS_H__

#include "utils_ch_implicit.h"

void set_boundary_ghost_nodes( AppCtx* user , PetscScalar*** uarray , PetscInt Mx , PetscInt My , PetscInt Mz , PetscInt i , PetscInt j , PetscInt k );

PetscReal reset_boundary_residual_values_for_neumann_bc( PetscReal*** uarray ,
							 PetscReal rhs_ijk ,
							 PetscReal udot_ijk ,
							 PetscInt Mx , PetscInt My , PetscInt Mz ,
							 PetscInt i , PetscInt j , PetscInt k );

PetscReal reset_boundary_residual_values_for_dirichlet_bottom_neumann_remainder_bc( PetscReal** uarray ,
										    PetscReal rhs_ij ,
										    PetscReal udot_ij ,
										    PetscInt Mx , PetscInt My ,
										    PetscInt i , PetscInt j );

PetscReal reset_boundary_residual_values_for_dirichlet_topandbottom_neumann_remainder_bc( PetscReal** uarray ,
											  PetscReal rhs_ij ,
											  PetscReal udot_ij ,
											  PetscInt Mx , PetscInt My ,
											  PetscInt i , PetscInt j );

#endif
