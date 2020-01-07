
static char help[] = "JFNK implicit solver for 2D CH with PETSc \n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscvec.h>
#include <petscviewer.h>
#include <petscsys.h>
#include <petscksp.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "utils_ch_implicit.h"
#include "boundary_conditions.h"
#include "initial_conditions.h"
#include "rhs_implicit.h"
#include "petsc_event_handling.h"
#include "temperature_dependence.h"

int main(int argc,char **argv) {
  
  TS             ts;                   /* nonlinear solver */
  Vec            u,r;                  /* solution, residual vectors */
  Mat            J,Jmf = NULL;         /* Jacobian matrices */
  PetscErrorCode ierr;
  DM             da;
  PetscReal      dt;
  SNES           snes;
  KSP            ksp;

  ierr = PetscInitialize(&argc,&argv,argv[1],help);if (ierr) return ierr;

  AppCtx         user = parse_petsc_options();

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create distributed array (DMDA) to manage parallel grid and vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */  
  DMDACreate3d( PETSC_COMM_WORLD, 
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,    // type of boundary nodes
                DMDA_STENCIL_BOX,                                                 // type of stencil
                11, 11, 11,                                                       // global dimns of array
                PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,                         // #procs in each dimn
                1,                                                                // DOF per node
                2,                                                                // Stencil width
                NULL, NULL, NULL, &da );
  DMSetFromOptions(da);
  DMSetUp(da);
  user.da = da;

  // Rescale value of L_omega to match the user-specified domain size
  PetscScalar L_domain = powf( user.Lx * user.Ly * user.Lz , 1./3. );
  user.L_omega        *= L_domain;

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Extract global vectors from DMDA;
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMCreateGlobalVector(da,&u);
  VecDuplicate(u,&r);
  VecDuplicate(u,&user.eps_2);
  VecDuplicate(u,&user.sigma);
  VecDuplicate(u,&user.temperature);
  VecDuplicate(u,&user.X);
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create timestepping solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  TSCreate(PETSC_COMM_WORLD,&ts);
  TSSetProblemType(ts,TS_NONLINEAR);
  TSSetDM(ts,da);
  TSSetIFunction(ts,r,FormIFunctionTEST,&user);
  TSSetMaxTime(ts,user.t_final);
  TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set initial conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  FormInitialSolution( u , user.temperature , &user );
  compute_eps2_and_sigma_from_temperature( &user );
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
  SNESSetFromOptions(snes);
  SNESGetKSP(snes,&ksp);
  KSPSetFromOptions(ksp);
  PetscOptionsView( NULL , PETSC_VIEWER_STDOUT_WORLD );

  /* Set directions and terminate flags for the two events */
  PetscInt       direction[2];
  PetscBool      terminate[2];
  direction[0] = 1; direction[1] = 1;
  terminate[0] = PETSC_FALSE; terminate[1] = PETSC_FALSE;
  TSSetEventHandler( ts , 2 , direction , terminate , EventFunction , PostEventFunction_ResetTemperatureGaussianProfile , (void*)&user );
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const std::string initial_soln = "c_" + std::to_string( 0.0 ).substr(0,6) + ".bin";
  PetscPrintf( PETSC_COMM_WORLD , "Logging initial solution at t = 0 seconds\n" );
  log_solution( u , initial_soln );

  TSSolve(ts,u);

  const std::string final_soln = "c_" + std::to_string( user.t_final ).substr(0,6) + ".bin";
  PetscPrintf( PETSC_COMM_WORLD , "Logging final solution at t = %5.4f seconds\n" , (double)user.t_final );
  log_solution( u , final_soln );
  
  PetscPrintf( PETSC_COMM_WORLD , "SIMULATION DONE\n\n" );
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  MatDestroy(&J);
  MatDestroy(&Jmf);
  VecDestroy(&u);
  VecDestroy(&r);
  VecDestroy(&user.eps_2);
  VecDestroy(&user.sigma);
  VecDestroy(&user.temperature);
  VecDestroy(&user.X);
  TSDestroy(&ts);
  DMDestroy(&da);

  PetscFinalize();
  return ierr;
}
