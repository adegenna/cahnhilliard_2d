
static char help[] = "Explicit solver for 2D CH with PETSc \n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscvec.h>
#include <petscviewer.h>
#include <petscsys.h>
#include <petscdmcomposite.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "utils_ch_explicit.h"
#include "boundary_conditions.h"
#include "initial_conditions.h"
#include "rhs_explicit.h"
#include "petsc_event_handling.h"
#include "temperature_dependence.h"

int main(int argc,char **argv) {

  PetscInitialize(&argc,&argv,argv[1],help);

  TS             ts;                   /* nonlinear solver */
  Vec            u,c,T,r,r_c,r_T;      /* solution, residual vectors */
  Vec            U_c , U_T;
  Mat            J,Jmf = NULL;         /* Jacobian matrices */
  PetscErrorCode ierr;
  DM             da_c , da_T , pack;
  PetscReal      dt;
  SNES           snes;
  KSP            ksp;

  AppCtx         user = parse_petsc_options();

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create distributed array (DMDA) to manage parallel grid and vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */  

  // DM for concentration c
  DMDACreate2d(PETSC_COMM_WORLD, 
               DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,    // type of boundary nodes
               DMDA_STENCIL_BOX,                // type of stencil
               11,11,                           // global dimns of array (will be overwritten from user-options)
               PETSC_DECIDE,PETSC_DECIDE,       // #procs in each dimn
               1,                               // DOF per node
               2,                               // Stencil width
               NULL,NULL,&da_c);
  DMSetFromOptions(da_c);
  DMSetOptionsPrefix(da_c,"c_");
  DMSetUp(da_c);

  // Get process ownership ranges so that you can link different physics with the same indices
  const PetscInt *lxc , *lyc;
  PetscInt *lxT , *lyT;
  PetscInt sizes_x , sizes_y;
  PetscInt nx , ny;
  DMDAGetOwnershipRanges( da_c , &lxc , &lyc , 0 );
  DMDAGetInfo( da_c , NULL, &nx,&ny,NULL, &sizes_x,&sizes_y,NULL, NULL,NULL,NULL,NULL,NULL,NULL );
  PetscMalloc1( sizes_x , &lxT );
  PetscMalloc1( sizes_y , &lyT );
  PetscMemcpy( lxT , lxc , sizes_x*sizeof(*lxc) );
  PetscMemcpy( lyT , lyc , sizes_y*sizeof(*lyc) );

  // DM for temperature T
  DMDACreate2d(PETSC_COMM_WORLD, 
               DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,    // type of boundary nodes
               DMDA_STENCIL_STAR,               // type of stencil
               nx,ny,                           // global dimns of array
               sizes_x,sizes_y,                 // #procs in each dimn
               1,                               // DOF per node
               2,                               // Stencil width
               lxT,lyT,&da_T);
  DMSetFromOptions(da_T);
  DMSetOptionsPrefix(da_T,"T_");
  DMSetUp(da_T);
  
  PetscFree(lxT);
  PetscFree(lyT);

  // Pack the DMs together into a single composite manager
  DMCompositeCreate( PETSC_COMM_WORLD , &pack );
  DMSetOptionsPrefix( pack , "pack_" );
  DMCompositeAddDM( pack , da_c );
  DMCompositeAddDM( pack , da_T );
  DMDASetFieldName( da_c , 0 , "c" );
  DMDASetFieldName( da_T , 0 , "T" );
  DMSetFromOptions( pack );
  
  // Rescale value of L_omega to match the user-specified domain size
  PetscScalar L_domain = sqrtf( user.Lx * user.Ly );
  user.L_omega        *= L_domain;

  user.da_c = da_c;
  user.da_T = da_T;
  user.pack = pack;
  
  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Extract global vectors from DMDA;
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMCreateGlobalVector( da_c , &c );
  DMCreateGlobalVector( da_T , &T );
  DMCreateGlobalVector( pack , &u );
  VecDuplicate( c , &r_c ); // Residual used for CH-only calculations
  VecDuplicate( T , &r_T ); // Residual used for thermal-only calculations
  VecDuplicate( u , &r );   // Residual used for coupled calculations
  VecDuplicate(c,&user.eps_2);
  VecDuplicate(c,&user.sigma);
  VecDuplicate(T,&user.temperature_source);
  VecDuplicate(c,&user.X);
  VecDuplicate(T,&user.dirichlet_bc_ch_array);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create timestepping solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  TSCreate(PETSC_COMM_WORLD,&ts);
  TSSetProblemType(ts,TS_NONLINEAR);
  TSSetMaxTime(ts,user.t_final);
  TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);
  TSSetTimeStep(ts,user.dt);

  // Initial solution
  FormInitialSolution( u , &user );
  
  // Thermal CH parameters
  compute_eps2_and_sigma_from_temperature( &user , u );
  
  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set type of physics
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  // Unpack the stuff you need
  DMCompositeGetAccess( pack , u , &U_c , &U_T );

  PetscErrorCode (*rhsFunctionExplicit)( TS ts , PetscReal t , Vec U , Vec F , void *ctx );
  DM  da_user;
  Vec U_user , r_user;
  
  if (user.physics.compare("ch") == 0) {
    // CH only
    
    da_user = da_c;
    U_user  = U_c;
    r_user  = r_c;
    rhsFunctionExplicit = FormRHS_CH;
    
  }
  
  else if (user.physics.compare("thermal") == 0) {
    // Thermal only

    da_user = da_T;
    U_user  = U_T;
    r_user  = r_T;
    rhsFunctionExplicit = FormRHS_thermal;
    
  }
  
  else if (user.physics.compare("coupled_ch_thermal") == 0) {
    // Coupled CH + thermal

    da_user = pack;
    U_user  = u;
    r_user  = r;
    rhsFunctionExplicit = FormRHS_CH_coupled;
    
  }

  else {
    // Incorrectly specified physics option
    
    PetscPrintf( PETSC_COMM_WORLD , "Error: physics option specified incorrectly ...\n\n" );

    return(0);

  }
  
  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set boundary condition function
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    // CH equation: c
  if ( user.boundary_ch.compare("neumann") == 0 ) // Neumann
    user.residualFunction_ch = reset_boundary_values_for_neumann_bc;

  else if ( user.boundary_ch.compare("dirichlet") == 0 ) // Dirichlet
    user.residualFunction_ch = reset_boundary_values_for_dirichlet_bc;
  
  else {
    // Incorrectly specified bc option
    
    PetscPrintf( PETSC_COMM_WORLD , "Error: boundary option specified incorrectly ...\n\n" );

    return(0);
    
  }
  
  // Thermal equation
  if ( user.boundary_thermal.compare("neumann") == 0 ) // Neumann
    user.residualFunction_thermal = reset_boundary_values_for_neumann_bc;

  else if ( user.boundary_thermal.compare("dirichlet") == 0 ) // Dirichlet
    user.residualFunction_thermal = reset_boundary_values_for_dirichlet_bc;
  
  else {
    // Incorrectly specified bc option
    
    PetscPrintf( PETSC_COMM_WORLD , "Error: boundary option specified incorrectly ...\n\n" );

    return(0);
    
  }
  
  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set time-stepping scheme
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  TSSetDM( ts , da_user );
  TSSetSolution( ts , U_user );  
  TSSetRHSFunction( ts , r_user , rhsFunctionExplicit , &user );

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set user options and event handling
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
  // User-options
  TSSetFromOptions(ts);
  SNESSetFromOptions(snes);
  PetscOptionsView( NULL , PETSC_VIEWER_STDOUT_WORLD );  
  
  // Setup event handling
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
  log_solution( U_c , initial_soln );

  TSSolve( ts , U_user );
    
  const std::string final_soln = "c_" + std::to_string( user.t_final ).substr(0,6) + ".bin";
  PetscPrintf( PETSC_COMM_WORLD , "Logging final solution at t = %5.4f seconds\n" , (double)user.t_final );
  log_solution( U_c , final_soln );
  
  PetscPrintf( PETSC_COMM_WORLD , "SIMULATION DONE\n\n" );
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMCompositeRestoreAccess( pack , u , &U_c , &U_T );
  MatDestroy(&J);
  MatDestroy(&Jmf);
  VecDestroy(&u);
  VecDestroy(&c);
  VecDestroy(&r_c);
  VecDestroy(&r_T);
  VecDestroy(&r);
  VecDestroy(&user.eps_2);
  VecDestroy(&user.sigma);
  VecDestroy(&user.temperature_source);
  VecDestroy(&user.X);
  VecDestroy(&U_user);
  VecDestroy(&r_user);
  TSDestroy(&ts);
  DMDestroy(&da_c);
  DMDestroy(&da_T);
  DMDestroy(&da_user);
  
  PetscFinalize();
  return ierr;
}
