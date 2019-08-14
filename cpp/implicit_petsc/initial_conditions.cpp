#include "initial_conditions.h"

PetscErrorCode FormInitialSolution(Vec U,void *ptr)
{
  AppCtx         *user=(AppCtx*)ptr;
  DM             da   =user->da;
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym,Mx,My;
  PetscScalar    **u;
  PetscReal      hx,hy,x,y,r;
  PetscRandom rng;
  PetscReal value_rng;

  PetscFunctionBeginUser;
  DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);

  hx = (user->Lx)/(PetscReal)(Mx-1);
  hy = (user->Ly)/(PetscReal)(My-1);

  /* Get pointers to vector data */
  DMDAVecGetArray(da,U,&u);

  /* Get local grid boundaries */
  DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);

  /* Compute function over the locally owned part of the grid */
  PetscRandomCreate(PETSC_COMM_WORLD,&rng);
  PetscRandomSetType(rng,PETSCRAND48);

  // Interior
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      PetscRandomGetValueReal(rng , &value_rng);
      u[j][i] = 0.005 * ( 2.0 * value_rng - 1.0 );
    }
  }
  
  PetscRandomDestroy(&rng);

  /* Restore vectors */
  DMDAVecRestoreArray(da,U,&u);
  PetscFunctionReturn(0);
}
