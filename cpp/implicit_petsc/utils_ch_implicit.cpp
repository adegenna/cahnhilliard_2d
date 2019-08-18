
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
  

  
  
  return user;
  
};
