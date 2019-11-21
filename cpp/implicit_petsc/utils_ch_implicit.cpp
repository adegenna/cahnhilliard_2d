#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <stdio.h>
#include "utils_ch_implicit.h"

AppCtx parse_petsc_options( ) {

  AppCtx         user;                 /* user-defined work context */
  
  /* Initialize user application context */

  user.da_c           = NULL;
  user.pack           = NULL;
  
  // Grid
  PetscOptionsGetReal(NULL,NULL,"-Lx",&user.Lx,NULL);
  PetscOptionsGetReal(NULL,NULL,"-Ly",&user.Ly,NULL);

  // Boundary conditions
  char tempfile_boundary[PETSC_MAX_PATH_LEN];
  PetscOptionsGetString(NULL,NULL,"-boundary",tempfile_boundary,sizeof(tempfile_boundary),NULL);
  user.boundary = std::string(tempfile_boundary);
  PetscOptionsGetReal(NULL,NULL,"-dirichlet_bc",&user.dirichlet_bc,NULL);
  
  // Temporal scheme
  PetscOptionsGetReal(NULL,NULL,"-t_final",&user.t_final,NULL);
  PetscOptionsGetReal(NULL,NULL,"-dt_check",&user.dt_check,NULL);
  PetscOptionsGetReal(NULL,NULL,"-dt_output",&user.dt_output,NULL);
  char tempfile_timestepper[PETSC_MAX_PATH_LEN];
  PetscOptionsGetString(NULL,NULL,"-time_stepper",tempfile_timestepper,sizeof(tempfile_timestepper),NULL);
  user.time_stepper = std::string(tempfile_timestepper);
  PetscOptionsGetReal(NULL,NULL,"-dt",&user.dt,NULL);

  // Thermal options
  char tempfile[PETSC_MAX_PATH_LEN];
  char tempfile_Tsource[PETSC_MAX_PATH_LEN];
  PetscOptionsGetString(NULL,NULL,"-initial_temperature_file",tempfile,sizeof(tempfile),NULL);
  user.initial_temperature_file = std::string(tempfile);
  PetscOptionsGetString(NULL,NULL,"-initial_temperature_source_file",tempfile_Tsource,sizeof(tempfile_Tsource),NULL);
  user.initial_temperature_source_file = std::string(tempfile_Tsource);
  PetscOptionsGetReal(NULL,NULL,"-D_T",&user.D_T,NULL);  

  // CH options
  PetscOptionsGetReal(NULL,NULL,"-CH_m",&user.m,NULL);
  char tempfile_U[PETSC_MAX_PATH_LEN];
  PetscOptionsGetString(NULL,NULL,"-initial_soln_file",tempfile_U,sizeof(tempfile_U),NULL);
  user.initial_soln_file = std::string(tempfile_U);
  char tempfile_physics[PETSC_MAX_PATH_LEN];
  PetscOptionsGetString(NULL,NULL,"-physics",tempfile_physics,sizeof(tempfile_physics),NULL);
  user.physics = std::string(tempfile_physics);

  return user;
  
};
