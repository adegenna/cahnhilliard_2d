#ifndef __BOUNDARY_CONDITIONS_H__
#define __BOUNDARY_CONDITIONS_H__

#include "utils_ch_explicit.h"

void set_boundary_ghost_nodes_dirichlet( AppCtx* user , PetscScalar** uarray , PetscInt Mx , PetscInt My , PetscInt i , PetscInt j );

void set_boundary_ghost_nodes_dirichlet_singleframe( AppCtx* user , PetscScalar** uarray , PetscInt Mx , PetscInt My , PetscInt i , PetscInt j );

void set_boundary_ghost_nodes_dirichlet_singleframe_thermal( AppCtx* user , PetscScalar** uarray , PetscInt Mx , PetscInt My , PetscInt i , PetscInt j );

void set_boundary_ghost_nodes_neumann( AppCtx* user , PetscScalar** uarray , PetscInt Mx , PetscInt My , PetscInt i , PetscInt j );

ThirteenPointStencil get_thirteen_point_stencil( AppCtx* user ,
						 PetscReal** uarray ,
						 PetscInt Mx , PetscInt My ,
						 PetscInt i , PetscInt j );

PetscReal reset_boundary_values_for_neumann_bc( PetscReal** uarray ,
						PetscReal** u_none ,
						PetscInt Mx , PetscInt My ,
						PetscInt i , PetscInt j );

PetscReal reset_boundary_values_for_dirichlet_bc( PetscReal** uarray ,
						  PetscReal** u_dirichlet ,
						  PetscInt Mx , PetscInt My ,
						  PetscInt i , PetscInt j );

PetscReal reset_boundary_residual_values_for_neumann_bc( PetscReal** uarray ,
							 PetscReal** u_none ,
							 PetscReal rhs_ij ,
							 PetscReal udot_ij ,
							 PetscInt Mx , PetscInt My ,
							 PetscInt i , PetscInt j );

PetscReal reset_boundary_residual_values_for_dirichlet_bc( PetscReal** uarray ,
							   PetscReal** u_dirichlet ,
							   PetscReal rhs_ij ,
							   PetscReal udot_ij ,
							   PetscInt Mx , PetscInt My ,
							   PetscInt i , PetscInt j );

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

PetscReal compute_residuals_no_explicit_boundary_resets(                                  PetscReal** uarray ,
											  PetscReal rhs_ij ,
											  PetscReal udot_ij ,
											  PetscInt Mx , PetscInt My ,
											  PetscInt i , PetscInt j );

ThirteenPointStencil apply_dirichlet_bc( AppCtx* user ,
					 PetscReal** uarray ,
					 PetscInt Mx , PetscInt My ,
					 PetscInt i , PetscInt j );

ThirteenPointStencil apply_neumann_bc( AppCtx* user ,
				       PetscReal** uarray ,
				       PetscInt Mx , PetscInt My ,
				       PetscInt i , PetscInt j );


#endif
