#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <stdio.h>
#include "utils_ch_implicit.h"

AppCtx parse_petsc_options( ) {

  AppCtx         user;                 /* user-defined work context */
  
  /* Initialize user application context */

  user.pack           = NULL;
  user.da_c           = NULL;
  user.da_phi         = NULL;
  
  // Grid
  PetscOptionsGetReal(NULL,NULL,"-Lx",&user.Lx,NULL);
  PetscOptionsGetReal(NULL,NULL,"-Ly",&user.Ly,NULL);
  PetscOptionsGetReal(NULL,NULL,"-Lz",&user.Lz,NULL);

  // Boundary conditions
  char tempfile_boundary[PETSC_MAX_PATH_LEN];
  PetscOptionsGetString(NULL,NULL,"-boundary",tempfile_boundary,sizeof(tempfile_boundary),NULL);
  user.boundary = std::string(tempfile_boundary);
  PetscOptionsGetReal(NULL,NULL,"-dirichlet_bc",&user.dirichlet_bc,NULL);
  
  // Temporal scheme
  PetscOptionsGetReal(NULL,NULL,"-t_final",&user.t_final,NULL);
  PetscOptionsGetReal(NULL,NULL,"-dt_check",&user.dt_check,NULL);
  PetscOptionsGetReal(NULL,NULL,"-dt_output",&user.dt_output,NULL);
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

DM createLinkedDA_starStencil3D( DM da_base , std::string fieldname ) {

  DM da_coupled;
  
  // Get process ownership ranges so that you can link different physics with the same indices
  const PetscInt *lxc , *lyc , *lzc;
  PetscInt sizes_x , sizes_y , sizes_z;
  PetscInt nx , ny , nz;
  DMDAGetOwnershipRanges( da_base , &lxc , &lyc , &lzc );
  DMDAGetInfo( da_base , NULL, &nx,&ny,&nz, &sizes_x,&sizes_y,&sizes_z, NULL,NULL,NULL,NULL,NULL,NULL );

  // Allocation for linking phi-physics
  PetscInt *lxPhi , *lyPhi , *lzPhi;
  PetscMalloc1( sizes_x , &lxPhi );
  PetscMalloc1( sizes_y , &lyPhi );
  PetscMalloc1( sizes_z , &lzPhi );
  PetscMemcpy( lxPhi , lxc , sizes_x*sizeof(*lxc) );
  PetscMemcpy( lyPhi , lyc , sizes_y*sizeof(*lyc) );
  PetscMemcpy( lzPhi , lzc , sizes_z*sizeof(*lzc) );

  // DM for phi
  DMDACreate3d( PETSC_COMM_WORLD, 
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, // type of boundary nodes
                DMDA_STENCIL_STAR,                                             // type of stencil
                nx,ny,nz,                                                      // global dimns of array
                sizes_x,sizes_y,sizes_z,                                       // #procs in each dimn
                1,                                                             // DOF per node
                1,                                                             // Stencil width
                lxPhi,lyPhi,lzPhi,&da_coupled);
  
  DMSetFromOptions(da_coupled);
  DMSetOptionsPrefix( da_coupled , fieldname.c_str() );
  DMSetUp(da_coupled);

  PetscFree(lxPhi);
  PetscFree(lyPhi);
  PetscFree(lzPhi);

  return da_coupled;

};
