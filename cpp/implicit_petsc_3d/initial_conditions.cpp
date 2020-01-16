#include "initial_conditions.h"
#include "temperature_dependence.h"
#include <stdio.h>
#include <petscviewerhdf5.h>
#include <petscdmcomposite.h>

PetscErrorCode FormInitialSolution(Vec U , void *ptr)
{
  AppCtx         *user = (AppCtx*)ptr;
  DM             pack  = user->pack;
  PetscErrorCode ierr;
  PetscInt       i,j,k,xs,ys,zs,xm,ym,zm;
  PetscScalar    ***u , ***phi;
  PetscReal      x,y,z,r;
  PetscRandom    rng;
  PetscReal      value_rng;
  Vec            U_c , U_phi , U_T;

  PetscFunctionBeginUser;

  // Unpack composite DM
  DMCompositeGetAccess( pack , U , &U_c , &U_phi , &U_T );
  
  // Interior
  PetscViewer viewer_U , viewer_T , viewer_Tsource , viewer_dirichlet_ch , viewer_dirichlet_thermal;
  PetscViewerBinaryOpen( PETSC_COMM_WORLD , user->initial_temperature_file.c_str()        , FILE_MODE_READ , &viewer_T );
  PetscViewerBinaryOpen( PETSC_COMM_WORLD , user->initial_temperature_source_file.c_str() , FILE_MODE_READ , &viewer_Tsource );
  PetscViewerBinaryOpen( PETSC_COMM_WORLD , user->initial_soln_file.c_str()               , FILE_MODE_READ , &viewer_U );
  VecLoad( U_T   , viewer_T );
  VecLoad( U_c   , viewer_U );
  VecLoad( user->temperature_source , viewer_Tsource );
  PetscViewerDestroy(&viewer_T);
  PetscViewerDestroy(&viewer_Tsource);
  PetscViewerDestroy(&viewer_U);

  if ( user->boundary_ch.compare("dirichlet") == 0 ) {
    PetscViewerBinaryOpen( PETSC_COMM_WORLD , user->dirichlet_ch_array_file.c_str()    , FILE_MODE_READ , &viewer_dirichlet_ch );
    VecLoad( user->dirichlet_bc_ch_array , viewer_dirichlet_ch );
    PetscViewerDestroy(&viewer_dirichlet_ch);
  }
  
  if ( user->boundary_thermal.compare("dirichlet") == 0 ) {
    PetscViewerBinaryOpen( PETSC_COMM_WORLD , user->dirichlet_thermal_array_file.c_str()    , FILE_MODE_READ , &viewer_dirichlet_thermal );
    VecLoad( user->dirichlet_bc_thermal_array , viewer_dirichlet_thermal );
    PetscViewerDestroy(&viewer_dirichlet_thermal);
  }

  // Populate phi
  DM da_c , da_phi , da_T;
  DMCompositeGetEntries( pack , &da_c , &da_phi , &da_T );
  DMDAVecGetArray( da_phi , U_phi , &phi );
  DMDAGetCorners( da_phi , &xs,&ys,&zs , &xm,&ym,&zm );
  for (k=zs; k<zs+zm; k++) {
    for (j=ys; j<ys+ym; j++) {
      for (i=xs; i<xs+xm; i++) {
        phi[k][j][i] = 1.0;
      }
    }
  }
  DMDAVecRestoreArray( da_phi , U_phi , &phi );

  // Repack everything
  DMCompositeRestoreAccess( pack , U , &U_c , &U_phi , &U_T );
  
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
