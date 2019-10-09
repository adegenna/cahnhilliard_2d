#include "petsc_event_handling.h"
#include <chrono>
#include <thread>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "utils_ch_implicit.h"
#include "temperature_dependence.h"
#include <petscvec.h>
#include <petscsys.h>
#include <petscis.h>
#include <petscviewer.h>
#include <petscviewerhdf5.h>
#include <petscdmcomposite.h>

void log_solution( Vec U , const std::string& outname ) {

  MPI_Comm       comm = PETSC_COMM_WORLD;
  PetscViewer    viewer;
  PetscObjectSetName((PetscObject) U, "state");

  // PetscViewerCreate( comm , &viewer );
  // PetscViewerSetType( viewer , PETSCVIEWERHDF5 );
  // PetscViewerFileSetMode( viewer , FILE_MODE_WRITE );
  
  // PetscViewerHDF5Open( comm , outname.c_str() , FILE_MODE_WRITE , &viewer );
  // PetscViewerHDF5PushGroup(viewer, "/stategroup");
  // VecView( U , viewer );
  // PetscViewerHDF5PopGroup(viewer);
  // PetscViewerDestroy( &viewer );

  PetscViewerBinaryOpen( comm , outname.c_str() , FILE_MODE_WRITE , &viewer );
  VecView( U , viewer );
  PetscViewerDestroy( &viewer );  

  return;

}

 
PetscErrorCode EventFunction( TS ts , PetscReal t , Vec U , PetscScalar *fvalue , void *ctx ) {

  AppCtx            *app=(AppCtx*)ctx;
  PetscErrorCode    ierr;

  // Event 1: stop for driver program 
  fvalue[0] = t - (app->dt_counter + 1) * app->dt_check;
  
  // Event 2: output solution
  fvalue[1] = t - (app->dt_output_counter + 1) * app->dt_output;

  return(0);

}

PetscErrorCode PostEventFunction_ResetTemperatureGaussianProfile(TS ts,PetscInt nevents,PetscInt event_list[],PetscReal t,Vec U,PetscBool forwardsolve,void* ctx) {

  AppCtx         *app = (AppCtx*)ctx;
  Vec            U_c , U_T;
  DM             da_c , da_T;
  
  DMCompositeGetEntries( app->pack , &da_c , &da_T );
  DMCompositeGetAccess( app->pack , U , &U_c , &U_T );

  for (int i=0; i<nevents; i++) {
    
    // Log solution
    if ( (event_list[i] == 1) && ( t < app->t_final ) ) {

      PetscPrintf( PETSC_COMM_WORLD , "Logging solution at t = %5.4f seconds\n" , (double)t );

      const std::string outname_C = "c_"           + std::to_string( (app->dt_output_counter + 1) * app->dt_output ).substr(0,6) + ".bin";
      const std::string outname_T = "temperature_" + std::to_string( (app->dt_output_counter + 1) * app->dt_output ).substr(0,6) + ".bin";
      
      log_solution( U_c , outname_C );
      log_solution( U_T , outname_T );
                    
      app->dt_output_counter += 1;

    }

    // Look for new input parameters from driver
    if ( (event_list[i] == 0) && ( t < app->t_final ) ) {

      const std::string name     = "T_" + std::to_string( (app->dt_counter + 1) * app->dt_check ).substr(0,6) + ".bin";
      const std::string petscout = "Attempting to read new T(amplitude,x_mean,y_mean,sigma) at t = %5.4f seconds from file " + name + "\n";
      PetscPrintf( PETSC_COMM_WORLD , petscout.c_str() , (double)t );

      // Output to file that you are ready for a new value
      const std::string outname  = "complete_" + std::to_string( (app->dt_counter + 1) * app->dt_check ).substr(0,6) + ".bin";
      MPI_Comm comm = PETSC_COMM_WORLD;
      PetscViewer viewer_out;
      PetscViewerCreate( comm , &viewer_out );
      PetscViewerBinaryOpen( comm , outname.c_str() , FILE_MODE_WRITE , &viewer_out );
      VecView( U , viewer_out );
      PetscViewerDestroy( &viewer_out );

      // PetscMPIInt myrank;
      // MPI_Comm_rank( comm , &myrank );
      // if ( myrank == 0 ) {
      // 	std::ofstream output( outname );  
      // 	output << "";
      // 	output.close();
      // }
      // PetscBarrier( NULL );
       
      // Wait until receive a new value of parameters from an input file
      while(true) {

	if (FILE *file = fopen(name.c_str(), "r") ) {

	  fclose(file);
        
	  std::this_thread::sleep_for( std::chrono::milliseconds(100) ); // Give the driver time to finish writing the new parameters

	  std::ifstream fin(name);
	  PetscScalar T_amp , T_x , T_y , T_sigma;
	  fin >> T_amp;
	  fin >> T_x;
	  fin >> T_y;
	  fin >> T_sigma;
        
	  compute_new_temperature_profile( app , U , T_amp , T_x , T_y , T_sigma );
          
	  app->dt_counter += 1;
        
	  PetscPrintf( PETSC_COMM_WORLD , "Changing (T_amp, T_x, T_y, T_sigma) at t = %5.4f seconds to ( %5.4f , %5.4f , %5.4f , %5.4f)\n" , (double)t , (double)T_amp , (double)T_x , (double)T_y , (double)T_sigma );

	  fin.close();        
	  break;
        
	}
	
      }

      // Recompute ch parameters based on new temperature
      compute_eps2_and_sigma_from_temperature( ctx , U );
    }
    
  }

  return(0);

}

void compute_new_temperature_profile( AppCtx* user , Vec U , PetscScalar T_amp , PetscScalar T_x , PetscScalar T_y , PetscScalar T_sigma  ) {

  DM             pack = user->pack;
  DM             da_c , da_T;
  PetscInt       i,j,xs,ys,xm,ym,Mx,My;
  PetscScalar    **T;
  PetscReal      x,y,r;
  Vec            U_c , U_T;

  PetscFunctionBeginUser;

  DMCompositeGetEntries( pack , &da_c , &da_T );
  
  DMDAGetInfo(da_T,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);

  DMDAVecGetArray( da_T , user->temperature_source , &T );

  DMDAGetCorners(da_T,&xs,&ys,NULL,&xm,&ym,NULL);

  // Interior
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      T[j][i]     = T_amp * PetscExpReal( -0.5 * ( (j-T_x)*(j-T_x) + (i-T_y)*(i-T_y) ) / (T_sigma * T_sigma) );
      T[j][i]     = std::max( T[j][i] , 0.5 );
    }
  }

  /* Restore vectors */
  DMDAVecRestoreArray( da_T , temperature_source , &T );

  return;


}
