
#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <stdio.h>
#include "utils_ch_implicit.h"

AppCtx parse_petsc_options( ) {

  AppCtx         user;                 /* user-defined work context */
  
  /* Initialize user application context */

  user.da           = NULL;

  // Grid
  PetscOptionsGetReal(NULL,NULL,"-Lx",&user.Lx,NULL);
  PetscOptionsGetReal(NULL,NULL,"-Ly",&user.Ly,NULL);

  // Boundary conditions
  PetscOptionsGetInt(NULL,NULL,"-boundary",&user.boundary,NULL);
  PetscOptionsGetReal(NULL,NULL,"-dirichlet_bc",&user.dirichlet_bc,NULL);
  
  // Temporal scheme
  PetscOptionsGetReal(NULL,NULL,"-t_final",&user.t_final,NULL);
  PetscOptionsGetReal(NULL,NULL,"-dt_check",&user.dt_check,NULL);
  PetscOptionsGetReal(NULL,NULL,"-dt_output",&user.dt_output,NULL);

  // Thermal options
  char tempfile[PETSC_MAX_PATH_LEN];
  PetscOptionsGetString(NULL,NULL,"-initial_temperature_file",tempfile,sizeof(tempfile),NULL);
  user.initial_temperature_file = std::string(tempfile);

  // CH options
  PetscOptionsGetReal(NULL,NULL,"-CH_m",&user.m,NULL);
  
  return user;
  
};
