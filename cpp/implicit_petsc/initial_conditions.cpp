#include "initial_conditions.h"
#include "temperature_dependence.h"
#include <stdio.h>
#include <iostream>
#include <fstream>

PetscErrorCode FormInitialSolution(Vec U , Vec Temperature , void *ptr)
{
  AppCtx         *user=(AppCtx*)ptr;
  DM             da   =user->da;
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym,Mx,My;
  PetscScalar    **u , **T;
  PetscReal      hx,hy,x,y,r;
  PetscRandom rng;
  PetscReal value_rng;

  PetscFunctionBeginUser;
  DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);

  hx = (user->Lx)/(PetscReal)(Mx-1);
  hy = (user->Ly)/(PetscReal)(My-1);

  /* Get pointers to vector data */
  DMDAVecGetArray(da,U,&u);
  DMDAVecGetArray(da,Temperature,&T);

  /* Get local grid boundaries */
  DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);

  /* Compute function over the locally owned part of the grid */
  PetscRandomCreate(PETSC_COMM_WORLD,&rng);
  PetscRandomSetType(rng,PETSCRAND48);

  // Interior
  std::ifstream Tin( user->initial_temperature_file );
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      PetscRandomGetValueReal(rng , &value_rng);
      u[j][i]     = 0.005 * ( 2.0 * value_rng - 1.0 );
      Tin >> T[j][i];
    } 
  }
  Tin.close();

  /* Restore vectors */
  DMDAVecRestoreArray(da,U,&u);
  DMDAVecRestoreArray(da,Temperature,&T);
  
  PetscRandomDestroy(&rng);


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
