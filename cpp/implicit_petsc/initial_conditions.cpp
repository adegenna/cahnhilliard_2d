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

  // Fill values
  PetscInt j_mir, i_mir;
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if ( (j >= 1) && (j <= My-2) && (i >= 1) && (i <= Mx-2) ) {
        // Interior
        PetscRandomGetValueReal(rng , &value_rng);
        u[j][i] = 0.005 * ( 2.0 * value_rng - 1.0 );
      }
      else {
        // Boundary
        if ( user->boundary == 0 )
          // Dirichlet
          u[j][i] = user->dirichlet_bc;
        else {
          // Neumann: just fill ghost cells with random values, and explicitly change residual at the end
          PetscRandomGetValueReal(rng , &value_rng);
          u[j][i] = 0.005 * ( 2.0 * value_rng - 1.0 );
        }
          
      }
      
    }
  }
  
  PetscRandomDestroy(&rng);

  /* Restore vectors */
  DMDAVecRestoreArray(da,U,&u);
  PetscFunctionReturn(0);
}
