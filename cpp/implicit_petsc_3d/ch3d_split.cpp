
static char help[] = "JFNK implicit solver for 2D CH with PETSc \n";

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
#include "utils_ch_implicit.h"
#include "boundary_conditions.h"
#include "initial_conditions.h"
#include "rhs_implicit.h"
#include "petsc_event_handling.h"
#include "temperature_dependence.h"

int main(int argc,char **argv) {

  PetscInitialize(&argc,&argv,argv[1],help);

  TS             ts;                          /* nonlinear solver */
  Vec            u,c,phi,T,r,r_c,r_phi,r_T;   /* solution, residual vectors */
  Vec            U_c , U_phi , U_T;
  Mat            J,Jmf = NULL;                /* Jacobian matrices */
  PetscErrorCode ierr;
  DM             da_c , da_phi , da_T , pack;
  PetscReal      dt;
  SNES           snes;
  KSP            ksp;

  AppCtx         user = parse_petsc_options();

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create distributed array (DMDA) to manage parallel grid and vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */  

  // DM for concentration c
  DMDACreate3d( PETSC_COMM_WORLD, 
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,   // type of boundary nodes
                DMDA_STENCIL_STAR,                                               // type of stencil
                1,1,1,                                                           // global dimns of array
                PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,                          // #procs in each dimn
                1,                                                               // DOF per node
                1,                                                               // Stencil width
                NULL, NULL, NULL, &da_c );
  DMSetFromOptions(da_c);
  DMSetOptionsPrefix(da_c,"c_");
  DMSetUp(da_c);

  da_phi = createLinkedDA_starStencil3D( da_c , "phi_" );
  da_T   = createLinkedDA_starStencil3D( da_phi , "T_" );

  // Pack the DMs together into a single composite manager
  DMCompositeCreate( PETSC_COMM_WORLD , &pack );
  DMSetOptionsPrefix( pack , "pack_" );
  DMCompositeAddDM( pack , da_c );
  DMCompositeAddDM( pack , da_phi );
  DMCompositeAddDM( pack , da_T );
  DMDASetFieldName( da_c , 0 , "c" );
  DMDASetFieldName( da_phi , 0 , "phi" );
  DMDASetFieldName( da_T , 0 , "T" );
  DMSetFromOptions( pack );
  
  // Rescale value of L_omega to match the user-specified domain size
  PetscScalar L_domain = powf( user.Lx * user.Ly * user.Lz , 1./3. );
  user.L_omega        *= L_domain;

  user.da_c   = da_c;
  user.da_phi = da_phi;
  user.da_T   = da_T;
  user.pack   = pack;
  
  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Extract global vectors from DMDA;
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMCreateGlobalVector( da_c   , &c );
  DMCreateGlobalVector( da_phi , &phi );
  DMCreateGlobalVector( da_T   , &T );
  DMCreateGlobalVector( pack   , &u );
  VecDuplicate( c   , &r_c ); // Residual used for CH-c calculations
  VecDuplicate( phi , &r_phi ); // Residual used for CH-phi calculations
  VecDuplicate( T , &r_T ); // Residual used for thermal-only calculations
  VecDuplicate( u , &r );   // Residual used for coupled calculations
  VecDuplicate( c   , &user.sigma );
  VecDuplicate( phi , &user.eps_2 );
  VecDuplicate( c   , &user.X );
  VecDuplicate( T   , &user.temperature_source );

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
  DMCompositeGetAccess( pack , u , &U_c , &U_phi , &U_T );

  PetscErrorCode (*rhsFunctionImplicit)( TS ts , PetscReal t , Vec U , Vec Udot , Vec F , void *ctx );
  DM  da_user;
  Vec U_user , r_user;
  
  if (user.physics.compare("ch") == 0) {
    // CH only
    
    da_user = pack;
    U_user  = u;
    r_user  = r;
    rhsFunctionImplicit = FormIFunction_CH_split;
    
  }

  else if (user.physics.compare("coupled_ch_thermal") == 0) {
    // Coupled CH-thermal
    
    da_user = pack;
    U_user  = u;
    r_user  = r;
    rhsFunctionImplicit = FormIFunction_CH_split_thermal;

  }

  else {
    // Incorrectly specified physics option
    
    PetscPrintf( PETSC_COMM_WORLD , "Error: physics option specified incorrectly ...\n\n" );

    return(0);

  }

  
  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set boundary condition function
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
  if ( user.boundary.compare("neumann") == 0 ) // Neumann
    user.residualFunction = reset_boundary_residual_values_for_neumann_bc;

  else if ( user.boundary.compare("dirichlet") == 0 ) // Dirichlet
    user.residualFunction = compute_residuals_no_explicit_boundary_resets;
  
  // else if ( user.boundary.compare("bottom_dirichlet_neumann_remainder") == 0 ) // Bottom dirichlet, rest Neumann
  //   user.residualFunction = reset_boundary_residual_values_for_dirichlet_bottom_neumann_remainder_bc;

  // else if ( user.boundary.compare("topandbottom_dirichlet_neumann_remainder") == 0 ) // Bottom/top dirichlet, rest Neumann
  //   user.residualFunction = reset_boundary_residual_values_for_dirichlet_topandbottom_neumann_remainder_bc;
      
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
  
  TSSetIFunction( ts , r_user , rhsFunctionImplicit , &user );
  DMSetMatType( da_user , MATAIJ );
  DMCreateMatrix( da_user , &J );
  TSGetSNES( ts , &snes );
  MatCreateSNESMF( snes , &Jmf );
  SNESSetJacobian( snes , Jmf , J , SNESComputeJacobianDefaultColor , 0 );

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set user options and event handling
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
  // User-options
  TSSetFromOptions(ts);
  SNESSetFromOptions(snes);
  SNESGetKSP(snes,&ksp);
  KSPSetFromOptions(ksp);
  PetscOptionsView( NULL , PETSC_VIEWER_STDOUT_WORLD );  
  
  // // Setup event handling
  // PetscInt       direction[3];
  // PetscBool      terminate[3];
  // direction[0] = 1;           direction[1] = 1;           direction[2] = 1;
  // terminate[0] = PETSC_FALSE; terminate[1] = PETSC_FALSE; terminate[2] = PETSC_FALSE;
  // TSSetEventHandler( ts , 3 , direction , terminate , EventFunction , PostEventFunction_ResetTemperatureGaussianProfile , (void*)&user );
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  const std::string initial_soln  = "c_" + std::to_string( 0.0 ).substr(0,6) + ".bin";
  PetscPrintf( PETSC_COMM_WORLD , "Logging initial CH solution at t = 0 seconds\n" );
  log_solution( U_c , initial_soln );
  const std::string initial_soln2 = "T_" + std::to_string( 0.0 ).substr(0,6) + ".bin";
  PetscPrintf( PETSC_COMM_WORLD , "Logging initial thermal solution at t = 0 seconds\n" );
  log_solution( U_T , initial_soln2 );

  TSSolve( ts , U_user );
    
  const std::string final_soln  = "c_" + std::to_string( user.t_final ).substr(0,6) + ".bin";
  PetscPrintf( PETSC_COMM_WORLD , "Logging final solution at t = %5.4f seconds\n" , (double)user.t_final );
  log_solution( U_c , final_soln );
  const std::string final_soln2 = "T_" + std::to_string( user.t_final ).substr(0,6) + ".bin";
  PetscPrintf( PETSC_COMM_WORLD , "Logging initial thermal solution at t = %5.4f seconds\n" , (double)user.t_final );
  log_solution( U_T , final_soln2 );

  PetscPrintf( PETSC_COMM_WORLD , "SIMULATION DONE\n\n" );
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMCompositeRestoreAccess( pack , u , &U_c , &U_phi );
  MatDestroy(&J);
  MatDestroy(&Jmf);
  VecDestroy(&u);
  VecDestroy(&c);
  VecDestroy(&phi);
  VecDestroy(&T);
  VecDestroy(&r_c);
  VecDestroy(&r_phi);
  VecDestroy(&r_T);
  VecDestroy(&r);
  VecDestroy(&user.eps_2);
  VecDestroy(&user.sigma);
  VecDestroy(&user.X);
  VecDestroy(&U_user);
  VecDestroy(&r_user);
  TSDestroy(&ts);
  DMDestroy(&da_c);
  DMDestroy(&da_phi);
  DMDestroy(&da_T);
  DMDestroy(&da_user);
  
  PetscFinalize();
  return ierr;
}
