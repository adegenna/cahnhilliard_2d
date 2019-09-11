
static char help[] = "JFNK implicit solver for 2D CH with PETSc \n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscvec.h>
#include <petscviewer.h>
#include <petscviewerhdf5.h>
#include <petscsys.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "utils_ch_implicit.h"

int main(int argc,char **argv) {
  
  Vec            U;                    /* solution, residual vectors */
  PetscErrorCode ierr;
  DM             da;
  SNES           snes;

  ierr = PetscInitialize(&argc,&argv,argv[1],help);if (ierr) return ierr;

  AppCtx         user = parse_petsc_options();

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create distributed array (DMDA) to manage parallel grid and vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */  
  DMDACreate2d(PETSC_COMM_WORLD, 
               DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,    // type of boundary nodes
               DMDA_STENCIL_BOX,                // type of stencil
               11,11,                           // global dimns of array
               PETSC_DECIDE,PETSC_DECIDE,       // #procs in each dimn
               1,                               // DOF per node
               2,                               // Stencil width
               NULL,NULL,&da);
  DMSetFromOptions(da);
  DMSetUp(da);
  user.da = da;

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Extract global vectors from DMDA;
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMCreateGlobalVector(da,&U);  

  MPI_Comm       comm = PETSC_COMM_WORLD;
  PetscInt num_solns  = (int) user.t_final/user.dt_output + 1;
  for (int i=0; i<num_solns; i++) {

    // Load original MPI solnfile
    const std::string solnfile = "c_" + std::to_string( user.dt_output * i ).substr(0,6) + ".bin";
    PetscViewer    viewer_input;
    PetscViewerBinaryOpen( comm , solnfile.c_str() , FILE_MODE_READ , &viewer_input );
    VecLoad( U , viewer_input );
    PetscViewerDestroy( &viewer_input );

    // Rewrite in serial to ASCII
    PetscViewer viewer_output;
    const std::string fileoutname = "c_" + std::to_string( user.dt_output * i ).substr(0,6) + ".ascii";
    PetscViewerCreate( comm , &viewer_output );
    PetscViewerSetType( viewer_output , PETSCVIEWERASCII );
    PetscViewerFileSetMode( viewer_output , FILE_MODE_WRITE );
    PetscViewerASCIIOpen( comm , fileoutname.c_str() , &viewer_output );
    VecView( U , viewer_output );
    PetscViewerDestroy( &viewer_output );

  }  

  VecDestroy(&U);
  DMDestroy(&da);

  PetscFinalize();
  return ierr;
}
