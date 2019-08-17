
static char help[] = "JFNK implicit solver for 2D CH with PETSc \n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include "utils_ch_implicit.h"
#include "boundary_conditions.h"
#include "initial_conditions.h"
#include "rhs_implicit.h"
#include "petsc_event_handling.h"

int main(int argc,char **argv) {
  
  TS             ts;                   /* nonlinear solver */
  Vec            u,r;                  /* solution, residual vectors */
  Vec            eps_2;
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
  VecDuplicate(u,&user.eps_2);

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
  FormInitialSolution( u , user.eps_2 , &user );
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
  PetscInt       direction[2];
  PetscBool      terminate[2];
  direction[0] = -1; direction[1] = 1;
  terminate[0] = PETSC_FALSE; terminate[1] = PETSC_FALSE;
  TSSetEventHandler( ts , 2 , direction , terminate , EventFunction , PostEventFunction , (void*)&user );
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  TSSolve(ts,u);
  PetscPrintf( PETSC_COMM_WORLD , "SIMULATION DONE\n\n" );
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
