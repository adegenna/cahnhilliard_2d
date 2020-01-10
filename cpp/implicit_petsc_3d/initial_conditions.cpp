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
  Vec            U_c , U_phi;

  PetscFunctionBeginUser;

  // Unpack composite DM
  DMCompositeGetAccess( pack , U , &U_c , &U_phi );
  
  // Interior
  PetscViewer viewer_U , viewer_T;
  //PetscViewerBinaryOpen( PETSC_COMM_WORLD , user->initial_temperature_file.c_str() , FILE_MODE_READ , &viewer_T );
  PetscViewerBinaryOpen( PETSC_COMM_WORLD , user->initial_soln_file.c_str()        , FILE_MODE_READ , &viewer_U );
  //VecLoad( Temperature , viewer_T );
  VecLoad( U_c           , viewer_U );
  //PetscViewerDestroy(&viewer_T);
  PetscViewerDestroy(&viewer_U);

  // Populate phi
  DM da_c , da_phi;
  DMCompositeGetEntries( pack , &da_c , &da_phi );
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
  DMCompositeRestoreAccess( pack , U , &U_c , &U_phi );
  
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
