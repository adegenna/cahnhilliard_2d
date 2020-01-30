#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <stdio.h>
#include "utils_ch_explicit.h"

AppCtx parse_petsc_options( ) {

  AppCtx         user;                 /* user-defined work context */
  
  /* Initialize user application context */

  user.da_c           = NULL;
  user.da_T           = NULL;
  user.pack           = NULL;
  
  // Grid
  PetscOptionsGetReal(NULL,NULL,"-Lx",&user.Lx,NULL);
  PetscOptionsGetReal(NULL,NULL,"-Ly",&user.Ly,NULL);

  // Boundary conditions
  char tempfile_boundary_ch[PETSC_MAX_PATH_LEN];
  PetscOptionsGetString(NULL,NULL,"-boundary_ch",tempfile_boundary_ch,sizeof(tempfile_boundary_ch),NULL);
  user.boundary_ch = std::string(tempfile_boundary_ch);
  char tempfile_boundary_thermal[PETSC_MAX_PATH_LEN];
  PetscOptionsGetString(NULL,NULL,"-boundary_thermal",tempfile_boundary_thermal,sizeof(tempfile_boundary_thermal),NULL);
  user.boundary_thermal = std::string(tempfile_boundary_thermal);
  PetscOptionsGetReal(NULL,NULL,"-dirichlet_bc",&user.dirichlet_bc,NULL);
  char tempfile_dirichlet_array[PETSC_MAX_PATH_LEN];
  PetscOptionsGetString(NULL,NULL,"-dirichlet_bc_thermal_array_file",tempfile_dirichlet_array,sizeof(tempfile_dirichlet_array),NULL);
  user.dirichlet_thermal_array_file = std::string(tempfile_dirichlet_array);
  char tempfile_dirichlet_array_ch[PETSC_MAX_PATH_LEN];
  PetscOptionsGetString(NULL,NULL,"-dirichlet_bc_ch_array_file",tempfile_dirichlet_array_ch,sizeof(tempfile_dirichlet_array_ch),NULL);
  user.dirichlet_ch_array_file = std::string(tempfile_dirichlet_array_ch);
  
  // Temporal scheme
  PetscOptionsGetReal(NULL,NULL,"-t_final",&user.t_final,NULL);
  PetscOptionsGetReal(NULL,NULL,"-dt_check",&user.dt_check,NULL);
  PetscOptionsGetReal(NULL,NULL,"-dt_output",&user.dt_output,NULL);
  PetscOptionsGetReal(NULL,NULL,"-dt",&user.dt,NULL);
  PetscOptionsGetReal(NULL,NULL,"-dt_thermal_reset",&user.dt_thermal_reset,NULL);

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

DM createLinkedDA_starStencil2D( DM da_base , std::string fieldname ) {

  DM da_coupled;
  
  // Get process ownership ranges so that you can link different physics with the same indices
  const PetscInt *lxc , *lyc;
  PetscInt sizes_x , sizes_y;
  PetscInt nx , ny;
  DMDAGetOwnershipRanges( da_base , &lxc , &lyc , 0 );
  DMDAGetInfo( da_base , NULL, &nx,&ny,NULL, &sizes_x,&sizes_y,NULL, NULL,NULL,NULL,NULL,NULL,NULL );

  // Allocation for linking phi-physics
  PetscInt *lxPhi , *lyPhi;
  PetscMalloc1( sizes_x , &lxPhi );
  PetscMalloc1( sizes_y , &lyPhi );
  PetscMemcpy( lxPhi , lxc , sizes_x*sizeof(*lxc) );
  PetscMemcpy( lyPhi , lyc , sizes_y*sizeof(*lyc) );
  
  // DM for phi
  DMDACreate2d(PETSC_COMM_WORLD, 
               DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,    // type of boundary nodes
               DMDA_STENCIL_STAR,               // type of stencil
               nx,ny,                           // global dimns of array (will be overwritten from user-options)
               sizes_x,sizes_y,                 // #procs in each dimn
               1,                               // DOF per node
               1,                               // Stencil width
               lxPhi,lyPhi,&da_coupled);
  
  DMSetFromOptions(da_coupled);
  DMSetOptionsPrefix( da_coupled , fieldname.c_str() );
  DMSetUp(da_coupled);

  PetscFree(lxPhi);
  PetscFree(lyPhi);

  return da_coupled;

};
