#include "initial_conditions.h"
#include "temperature_dependence.h"
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
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym;
  PetscScalar    **u , **phi , **T;
  PetscReal      x,y,r;
  PetscRandom    rng;
  PetscReal      value_rng;
  Vec            U_c , U_phi , U_T;

  PetscFunctionBeginUser;

  // Unpack composite DM
  DM da_c , da_phi , da_T;  
  if (user->physics.compare("coupled_ch_thermal") == 0) {
    DMCompositeGetEntries( pack , &da_c , &da_phi , &da_T );
    DMCompositeGetAccess(  pack , U     , &U_c    , &U_phi    , &U_T );
  }
  else if (user->physics.compare("ch") == 0) {
    DMCompositeGetEntries( pack , &da_c    , &da_phi );
    DMCompositeGetAccess(  pack , U , &U_c , &U_phi );
    da_T = user->da_T;
    DMGetGlobalVector( da_T , &U_T );
  }
  
  // Load initial concentration and temperature from file
  PetscViewer viewer_T , viewer_U , viewer_Tsource; 
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

  // Populate phi
  DMDAVecGetArray( da_phi , U_phi , &phi );
  DMDAGetCorners( da_phi , &xs,&ys , NULL , &xm,&ym , NULL);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      phi[j][i] = 1.0;
    }
  }
  DMDAVecRestoreArray( da_phi , U_phi , &phi );

  // Repack everything
  if (user->physics.compare("coupled_ch_thermal") == 0) {
    DMCompositeRestoreAccess( pack , U , &U_c , &U_phi , &U_T );
  }
  else if (user->physics.compare("ch") == 0) {
    DMCompositeRestoreAccess( pack , U , &U_c , &U_phi );
    DMRestoreGlobalVector( da_T , &U_T );
  }
  
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


/* Compute function over the locally owned part of the grid */
// DMDAVecGetArray(da,U,&u);
// DMDAVecGetArray(da,Temperature,&T);
// DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);
// PetscRandomCreate(PETSC_COMM_WORLD,&rng);
// PetscRandomSetType(rng,PETSCRAND48);

// for (j=ys; j<ys+ym; j++) {
//   for (i=xs; i<xs+xm; i++) {
//     PetscRandomGetValueReal(rng , &value_rng);
//     u[j][i]     = 0.005 * ( 2.0 * value_rng - 1.0 );
//     T[j][i]     = 1.0;
//   } 
// }
// DMDAVecRestoreArray(da,U,&u);
// DMDAVecRestoreArray(da,Temperature,&T);  
// PetscRandomDestroy(&rng);
