
static char help[] = "JFNK implicit solver for 2D CH with PETSc \n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <stdio.h>
#include "utils_ch_implicit.h"
#include "boundary_conditions.h"
#include "initial_conditions.h"
#include "rhs_implicit.h"

PetscErrorCode EventFunction( TS ts , PetscReal t , Vec U , PetscScalar *fvalue , void *ctx ) {

  AppCtx            *app=(AppCtx*)ctx;
  PetscErrorCode    ierr;
  
  /* Event for ball height */
  if ( (t - (app->dt_counter + 1) * app->dt_check) < 0)
    fvalue[0] = 1;
  else
    fvalue[0] = 0;
  
  return(0);

}

PetscErrorCode PostEventFunction(TS ts,PetscInt nevents,PetscInt event_list[],PetscReal t,Vec U,PetscBool forwardsolve,void* ctx) {

  AppCtx         *app=(AppCtx*)ctx;
  PetscReal       m;

  if ( event_list[0] == 0 ) {
    app->m += 0.1;
    app->dt_counter += 1;
    PetscPrintf( PETSC_COMM_WORLD , "Changing m at t = %5.2f seconds to %5.2f\n" , (double)t , (double)app->m );
  }

  return(0);

}


int main(int argc,char **argv) {
  
  TS             ts;                   /* nonlinear solver */
  Vec            u,r;                  /* solution, residual vectors */
  Mat            J,Jmf = NULL;         /* Jacobian matrices */
  PetscErrorCode ierr;
  DM             da;
  PetscReal      dt;
  SNES           snes;

  PetscInitialize( NULL , NULL , argv[1] , help );
  
  AppCtx         user = parse_petsc_options();

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create distributed array (DMDA) to manage parallel grid and vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMDACreate2d(PETSC_COMM_WORLD, 
		      DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,      // type of ghost nodes
 		      DMDA_STENCIL_BOX,                        // type of stencil
		      11,11,                                   // global dimns of array
		      PETSC_DECIDE,PETSC_DECIDE,               // #procs in each dimn
		      1,                                       // DOF per node
		      2,                                       // Stencil width
		      NULL,NULL,&da);
  DMSetFromOptions(da);
  DMSetUp(da);
  user.da = da;

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Extract global vectors from DMDA;
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMCreateGlobalVector(da,&u);
  VecDuplicate(u,&r);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create timestepping solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  TSCreate(PETSC_COMM_WORLD,&ts);
  TSSetProblemType(ts,TS_NONLINEAR);
  TSSetDM(ts,da);
  TSSetIFunction(ts,r,FormIFunction,&user);
  TSSetMaxTime(ts,user.t_final);
  TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set initial conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  FormInitialSolution(u,&user);
  TSSetSolution(ts,u);
  dt   = 0.005;
  TSSetTimeStep(ts,dt);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set Jacobian evaluation routine
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr  = DMSetMatType(da,MATAIJ);
  ierr  = DMCreateMatrix(da,&J);
  TSGetSNES(ts,&snes);
  MatCreateSNESMF(snes,&Jmf);
  SNESSetJacobian(snes,Jmf,J,SNESComputeJacobianDefaultColor,0);

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Sets various TS parameters from user options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  TSSetFromOptions(ts);
  PetscOptionsView( NULL , PETSC_VIEWER_STDOUT_WORLD );

  /* Set directions and terminate flags for the two events */
  PetscInt       direction[1];
  PetscBool      terminate[1];
  direction[0] = -1;
  terminate[0] = PETSC_FALSE;
  TSSetEventHandler( ts , 1 , direction , terminate , EventFunction , PostEventFunction , (void*)&user );
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  TSSolve(ts,u);
  //printf("Simulation done, press enter to continue...\n");
  //getchar();
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  MatDestroy(&J);
  MatDestroy(&Jmf);
  VecDestroy(&u);
  VecDestroy(&r);
  TSDestroy(&ts);
  DMDestroy(&da);

  PetscFinalize();
  return ierr;
}
