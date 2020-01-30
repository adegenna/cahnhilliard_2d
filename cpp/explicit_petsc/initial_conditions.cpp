#include "initial_conditions.h"
#include "temperature_dependence.h"
#include "boundary_conditions.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <petscviewerhdf5.h>
#include <petscdmcomposite.h>

PetscErrorCode FormInitialSolution(Vec U , void *ptr)
{

  // Inputs:
  // Vec U : global packed vector (i.e. a PETSc concatentation of c & T)
  // void *ptr: the AppCtx struct defining various things for the simulation
  
  AppCtx         *user = (AppCtx*)ptr;
  DM             pack  = user->pack; // DM for coupled-CH
  DM             da_c , da_T;
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym,Mx,My;
  PetscScalar    **u , **T;
  PetscReal      hx,hy,x,y,r;
  PetscRandom    rng;
  PetscReal      value_rng;
  Vec            U_c , U_T;

  PetscFunctionBeginUser;

  // Unpack composite DM
  DMCompositeGetEntries( pack , &da_c , &da_T );
  DMCompositeGetAccess( pack , U , &U_c , &U_T );
  
  DMDAGetInfo(da_c,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);

  // NOTE: these CH eqns are dimensionless with domain length scale = 1. Physical domain size shows up in L_omega.
  hx = 1.0/(PetscReal)(Mx-1);
  hy = 1.0/(PetscReal)(My-1);

  // Load initial concentration and temperature from file
  PetscViewer viewer_T , viewer_Tsource , viewer_U , viewer_dirichlet_ch , viewer_dirichlet_thermal;
  MPI_Comm comm = PETSC_COMM_WORLD;
  PetscViewer    viewer_out;
  PetscViewerCreate( comm , &viewer_out );
  PetscViewerBinaryOpen( PETSC_COMM_WORLD , user->initial_temperature_file.c_str() , FILE_MODE_READ , &viewer_T );
  PetscViewerBinaryOpen( PETSC_COMM_WORLD , user->initial_temperature_source_file.c_str() , FILE_MODE_READ , &viewer_Tsource );
  PetscViewerBinaryOpen( PETSC_COMM_WORLD , user->initial_soln_file.c_str()        , FILE_MODE_READ , &viewer_U );
  VecLoad( U_T , viewer_T );
  VecLoad( U_c , viewer_U );
  VecLoad( user->temperature_source , viewer_Tsource );
  PetscViewerDestroy(&viewer_T);
  PetscViewerDestroy(&viewer_Tsource);
  PetscViewerDestroy(&viewer_U);

  // Set boundary values
  Vec local_c, local_T , local_c_dirichlet, local_T_dirichlet;
  PetscScalar **u_dirichlet , **T_dirichlet;
  DMDAGetCorners(da_c,&xs,&ys,NULL,&xm,&ym,NULL);
  DMGetLocalVector( da_c , &local_c );
  DMGlobalToLocalBegin( da_c , U_c , INSERT_VALUES , local_c );
  DMGlobalToLocalEnd(   da_c , U_c , INSERT_VALUES , local_c );
  DMDAVecGetArray( da_c , local_c , &u );
  
  if ( user->boundary_ch.compare("dirichlet") == 0 ) {
    PetscViewerBinaryOpen( PETSC_COMM_WORLD , user->dirichlet_ch_array_file.c_str()    , FILE_MODE_READ , &viewer_dirichlet_ch );
    VecLoad( user->dirichlet_bc_ch_array , viewer_dirichlet_ch );
    PetscViewerDestroy(&viewer_dirichlet_ch);
    DMGetLocalVector(     da_c , &local_c_dirichlet );
    DMGlobalToLocalBegin( da_c , user->dirichlet_bc_ch_array , INSERT_VALUES , local_c_dirichlet );
    DMGlobalToLocalEnd(   da_c , user->dirichlet_bc_ch_array , INSERT_VALUES , local_c_dirichlet );
    DMDAVecGetArray(      da_c , local_c_dirichlet , &u_dirichlet );
    for (int j=ys; j<ys+ym; j++) {
      for (int i=xs; i<xs+xm; i++) {
        if ( i <= 1 || j <= 1 || i >= (Mx-2) || j >= (My-2) )
          u[j][i] = u_dirichlet[j][i];
      }
    }
  }

  else if ( user->boundary_ch.compare("neumann") == 0 ) {
    for (int j=ys; j<ys+ym; j++) {
      for (int i=xs; i<xs+xm; i++) {
        u[j][i] = reset_boundary_rhs_values_for_neumann_bc( u , Mx , My , i , j );
      } 
    }
  }

  DMDAVecRestoreArray(  da_c , local_c , &u );
  DMLocalToGlobalBegin( da_c , local_c , INSERT_VALUES , U_c );
  DMLocalToGlobalEnd(   da_c , local_c , INSERT_VALUES , U_c );
  DMRestoreLocalVector( da_c , &local_c );

  DMDAGetCorners(da_T,&xs,&ys,NULL,&xm,&ym,NULL);
  DMGetLocalVector(     da_T , &local_T );
  DMGlobalToLocalBegin( da_T , U_T , INSERT_VALUES , local_T );
  DMGlobalToLocalEnd(   da_T , U_T , INSERT_VALUES , local_T );
  DMDAVecGetArray(      da_T , local_T , &T );
  
  if ( user->boundary_thermal.compare("dirichlet") == 0 ) {
    PetscViewerBinaryOpen( PETSC_COMM_WORLD , user->dirichlet_thermal_array_file.c_str()    , FILE_MODE_READ , &viewer_dirichlet_thermal );
    VecLoad( user->dirichlet_bc_thermal_array , viewer_dirichlet_thermal );
    PetscViewerDestroy(&viewer_dirichlet_thermal);
    DMGetLocalVector(     da_T , &local_T_dirichlet );
    DMGlobalToLocalBegin( da_T , user->dirichlet_bc_thermal_array , INSERT_VALUES , local_T_dirichlet );
    DMGlobalToLocalEnd(   da_T , user->dirichlet_bc_thermal_array , INSERT_VALUES , local_T_dirichlet );
    DMDAVecGetArray(      da_T , local_T_dirichlet , &T_dirichlet );
    for (int j=ys; j<ys+ym; j++) {
      for (int i=xs; i<xs+xm; i++) {
        if ( i <= 1 || j <= 1 || i >= (Mx-2) || j >= (My-2) )
          T[j][i] = T_dirichlet[j][i];
      }
    }

  }
  
  else if ( user->boundary_thermal.compare("neumann") == 0 ) {
    for (int j=ys; j<ys+ym; j++) {
      for (int i=xs; i<xs+xm; i++) {
        T[j][i] = reset_boundary_rhs_values_for_neumann_bc( T , Mx , My , i , j );
      } 
    }
  }

  DMDAVecRestoreArray(  da_T , local_T , &u );
  DMLocalToGlobalBegin( da_T , local_T , INSERT_VALUES , U_T );
  DMLocalToGlobalEnd(   da_T , local_T , INSERT_VALUES , U_T );
  DMRestoreLocalVector( da_T , &local_T );
  
  DMCompositeRestoreAccess( pack , U , &U_c , &U_T );
  
  // Compute temperature-dependent polymer limiters
  user->eps2_min = compute_eps2_from_chparams( user->X_max ,
                                               user->L_kuhn ,
                                               user->m ,
                                               user->L_omega );

  user->eps2_max = compute_eps2_from_chparams( user->X_min ,
                                               user->L_kuhn ,
                                               user->m ,
                                               user->L_omega );

  user->sigma_min = compute_sigma_from_chparams( user->X_max ,
                                                 user->L_kuhn ,
                                                 user->m ,
                                                 user->L_omega ,
                                                 user->N );

  user->sigma_max = compute_sigma_from_chparams( user->X_min ,
                                                 user->L_kuhn ,
                                                 user->m ,
                                                 user->L_omega ,
                                                 user->N );
  
  PetscFunctionReturn(0);
}
