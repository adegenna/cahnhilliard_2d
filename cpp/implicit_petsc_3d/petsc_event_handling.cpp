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

PetscErrorCode log_solution( Vec U , const std::string& outname ) {

  PetscFunctionBeginUser;
  
  MPI_Comm       comm = PETSC_COMM_WORLD;
  PetscViewer    viewer;
  PetscObjectSetName((PetscObject) U, "state");

  PetscViewerBinaryOpen( comm , outname.c_str() , FILE_MODE_WRITE , &viewer );
  VecView( U , viewer );
  PetscViewerDestroy( &viewer );  
  
  PetscFunctionReturn(0);

}

 
PetscErrorCode EventFunction( TS ts , PetscReal t , Vec U , PetscScalar *fvalue , void *ctx ) {

  AppCtx            *app=(AppCtx*)ctx;
  PetscErrorCode    ierr;
  
  PetscFunctionBeginUser;

  // Event 1: stop for driver program 
  fvalue[0] = t - (app->dt_counter + 1) * app->dt_check;
  
  // Event 2: output solution
  fvalue[1] = t - (app->dt_output_counter + 1) * app->dt_output;
  
  // Event 3: recalculate thermal-dependent parameters each new solution timestep
  fvalue[2] = t - (app->dt_thermal_counter + 1) * app->dt_thermal_reset;

  PetscFunctionReturn(0);

}

PetscErrorCode PostEventFunction_ResetM(TS ts,PetscInt nevents,PetscInt event_list[],PetscReal t,Vec U,PetscBool forwardsolve,void* ctx) {

  AppCtx         *app=(AppCtx*)ctx;
  PetscReal       m, m_new;

  PetscFunctionBeginUser;
  
  for (int i=0; i<nevents; i++) {
    
    // Log solution
    if ( (event_list[i] == 1) && ( t < app->t_final ) ) {

      PetscPrintf( PETSC_COMM_WORLD , "Logging solution at t = %5.4f seconds\n" , (double)t );

      const std::string outname = "c_" + std::to_string( (app->dt_output_counter + 1) * app->dt_output ).substr(0,6) + ".bin";

      log_solution( U , outname );
      
      app->dt_output_counter += 1;

    }

    // Look for new input parameters from driver
    if ( (event_list[i] == 0) && ( t < app->t_final ) ) {

      const std::string name     = "m_" + std::to_string( (app->dt_counter + 1) * app->dt_check ).substr(0,6) + ".out";
      const std::string petscout = "Attempting to read new m at t = %5.4f seconds from file " + name + "\n";
      PetscPrintf( PETSC_COMM_WORLD , petscout.c_str() , (double)t );

      // Output to file that you are ready for a new value
      const std::string outname  = "complete_" + std::to_string( (app->dt_counter + 1) * app->dt_check ).substr(0,6) + ".out";
      std::ofstream fout(outname);
      fout << "complete\n";
      fout.close();

      // Wait until receive a new value of m from an input file
      while(true) {

        if (FILE *file = fopen(name.c_str(), "r") ) {

          fclose(file);
        
          std::this_thread::sleep_for( std::chrono::milliseconds(100) ); // Give the driver time to finish writing the new parameters

          std::ifstream fin(name);
          PetscReal m_new;
          fin >> m_new;
        
          app->m           = m_new;
          app->dt_counter += 1;
        
          PetscPrintf( PETSC_COMM_WORLD , "Changing m at t = %5.4f seconds to m = %5.4f\n" , (double)t , (double)app->m );

          fin.close();        
          break;
        
        }

      }
      
      // Recompute ch parameters based on new temperature
      //compute_eps2_and_sigma_from_temperature( ctx );
    }
    
  }

  PetscFunctionReturn(0);

}

PetscErrorCode PostEventFunction_RecomputeThermalProperties(TS ts,PetscInt nevents,PetscInt event_list[],PetscReal t,Vec U,PetscBool forwardsolve,void* ctx) {

  AppCtx         *app = (AppCtx*)ctx;
  
  for (int i=0; i<nevents; i++) {

    if ( (event_list[i] == 2) && ( t < app->t_final ) ) {

      PetscPrintf( PETSC_COMM_WORLD , "Recalculating thermal properties t = %5.4f seconds\n" , (double)t );
      
      compute_eps2_and_sigma_from_temperature( ctx , U );

      app->dt_thermal_counter += 1;

    }

    // Log solution
    if ( (event_list[i] == 1) && ( t < app->t_final ) ) {

      PetscPrintf( PETSC_COMM_WORLD , "Logging solution at t = %5.4f seconds\n" , (double)t );

      const std::string outname = "c_" + std::to_string( (app->dt_output_counter + 1) * app->dt_output ).substr(0,6) + ".bin";

      log_solution( U , outname );
      
      app->dt_output_counter += 1;

    }

    
  }
  
  PetscFunctionReturn(0);

}

PetscErrorCode PostEventFunction_ResetTemperatureGaussianProfile(TS ts,PetscInt nevents,PetscInt event_list[],PetscReal t,Vec U,PetscBool forwardsolve,void* ctx) {

  AppCtx         *app = (AppCtx*)ctx;

  PetscFunctionBeginUser;

  for (int i=0; i<nevents; i++) {
    
    // Log solution
    if ( (event_list[i] == 1) && ( t < app->t_final ) ) {

      PetscPrintf( PETSC_COMM_WORLD , "Logging solution at t = %5.4f seconds\n" , (double)t );

      const std::string outname = "c_" + std::to_string( (app->dt_output_counter + 1) * app->dt_output ).substr(0,6) + ".bin";
      
      log_solution( U , outname );

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
       
      // Wait until receive a new value of parameters from an input file
      while(true) {

	if (FILE *file = fopen(name.c_str(), "r") ) {

	  fclose(file);
        
	  std::this_thread::sleep_for( std::chrono::milliseconds(100) ); // Give the driver time to finish writing the new parameters

	  std::ifstream fin(name);
	  PetscScalar T_amp , T_x , T_y , T_z , T_sigma;
	  fin >> T_amp;
	  fin >> T_x;
	  fin >> T_y;
	  fin >> T_z;
	  fin >> T_sigma;
        
	  compute_new_temperature_profile( app , T_amp , T_x , T_y , T_z , T_sigma );
          
	  app->dt_counter += 1;
        
	  PetscPrintf( PETSC_COMM_WORLD , "Changing (T_amp, T_x, T_y, T_z, T_sigma) at t = %5.4f seconds to ( %5.4f , %5.4f , %5.4f , %5.4f , %5.4f)\n" , (double)t , (double)T_amp , (double)T_x , (double)T_y , (double)T_z , (double)T_sigma );

	  fin.close();        
	  break;
        
	}
	
      }

      // Recompute ch parameters based on new temperature
      //compute_eps2_and_sigma_from_temperature( ctx );
    }
    
  }
  
  PetscFunctionReturn(0);

}

PetscErrorCode compute_new_temperature_profile( AppCtx* user , PetscScalar T_amp , PetscScalar T_x , PetscScalar T_y , PetscScalar T_z , PetscScalar T_sigma  ) {

  DM             pack   =user->pack;
  DM             da_c , da_phi , da_T;
  PetscInt       i,j,k,xs,ys,zs,xm,ym,zm,Mx,My,Mz;
  PetscScalar    ***T;
  PetscReal      x,y,z,r;
  Vec            U , U_c , U_phi , U_T;
  
  PetscFunctionBeginUser;

  DMCompositeGetEntries( pack , &da_c , &da_phi , &da_T );
  DMCompositeGetAccess(  pack , U     , &U_c    , &U_phi , &U_T );

  DMDAGetInfo( da_c ,
	       PETSC_IGNORE ,
	       &Mx , &My , &Mz ,
	       PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);

  /* Get pointers to vector data */
  DMDAVecGetArray(da_T,U_T,&T);

  /* Get local grid boundaries */
  DMDAGetCorners( da_c ,
		  &xs , &ys , &zs ,
		  &xm , &ym , &zm );

  // Interior
  for (k=zs; k<zs+zm; k++) {
    for (j=ys; j<ys+ym; j++) {
      for (i=xs; i<xs+xm; i++) {
	T[k][j][i]     = T_amp * PetscExpReal( -0.5 * ( (k-T_z)*(k-T_z) + (j-T_x)*(j-T_x) + (i-T_y)*(i-T_y) ) / (T_sigma * T_sigma) );
	T[k][j][i]     = std::max( T[k][j][i] , 0.5 );
      }
    }
  }
  
  /* Restore vectors */
  DMDAVecRestoreArray(da_T,U_T,&T);
  
  PetscFunctionReturn(0);


}
