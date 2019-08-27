#include "petsc_event_handling.h"
#include <chrono>
#include <thread>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "utils_ch_implicit.h"
#include "temperature_dependence.h"

PetscErrorCode EventFunction( TS ts , PetscReal t , Vec U , PetscScalar *fvalue , void *ctx ) {

  AppCtx            *app=(AppCtx*)ctx;
  PetscErrorCode    ierr;

  // Event 1: stop for driver program 
  fvalue[0] = t - (app->dt_counter + 1) * app->dt_check;
  
  // Event 2: output solution
  fvalue[1] = t - (app->dt_output_counter + 1) * app->dt_output;

  return(0);

}

PetscErrorCode PostEventFunction_ResetM(TS ts,PetscInt nevents,PetscInt event_list[],PetscReal t,Vec U,PetscBool forwardsolve,void* ctx) {

  AppCtx         *app=(AppCtx*)ctx;
  PetscReal       m, m_new;

  for (int i=0; i<nevents; i++) {
    
    // Log solution
    if ( (event_list[i] == 1) && ( t < app->t_final ) ) {

      PetscPrintf( PETSC_COMM_WORLD , "Logging solution at t = %5.4f seconds\n" , (double)t );

      const std::string outname = "c_" + std::to_string( (app->dt_output_counter + 1) * app->dt_output ).substr(0,6) + ".dat";

      const PetscScalar *u;

      PetscViewer    viewer;
      PetscViewerCreate( PETSC_COMM_WORLD , &viewer );
      PetscViewerSetType( viewer , PETSCVIEWERASCII );
      PetscViewerFileSetMode( viewer , FILE_MODE_WRITE );

      PetscViewerASCIIOpen( PETSC_COMM_WORLD , outname.c_str() , &viewer );
      VecView( U , viewer );
    
      PetscViewerDestroy( &viewer );
    
      app->dt_output_counter += 1;

    }

    // Look for new input parameters from driver
    if ( (event_list[i] == 0) && ( t < app->t_final ) ) {

      const std::string name     = "m_" + std::to_string( (app->dt_counter + 1) * app->dt_check ).substr(0,6) + ".dat";
      const std::string petscout = "Attepting to read new m at t = %5.4f seconds from file " + name + "\n";
      PetscPrintf( PETSC_COMM_WORLD , petscout.c_str() , (double)t );

      // Output to file that you are ready for a new value
      const std::string outname  = "complete_" + std::to_string( (app->dt_counter + 1) * app->dt_check ).substr(0,6) + ".dat";
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
      compute_eps2_and_sigma_from_temperature( ctx );
    }
    
  }

  return(0);

}

PetscErrorCode PostEventFunction_ResetTemperatureGaussianProfile(TS ts,PetscInt nevents,PetscInt event_list[],PetscReal t,Vec U,PetscBool forwardsolve,void* ctx) {

  AppCtx         *app=(AppCtx*)ctx;
  PetscReal       m, m_new;

  for (int i=0; i<nevents; i++) {
    
    // Log solution
    if ( (event_list[i] == 1) && ( t < app->t_final ) ) {

      PetscPrintf( PETSC_COMM_WORLD , "Logging solution at t = %5.4f seconds\n" , (double)t );

      const std::string outname = "c_" + std::to_string( (app->dt_output_counter + 1) * app->dt_output ).substr(0,6) + ".dat";

      const PetscScalar *u;

      PetscViewer    viewer;
      PetscViewerCreate( PETSC_COMM_WORLD , &viewer );
      PetscViewerSetType( viewer , PETSCVIEWERASCII );
      PetscViewerFileSetMode( viewer , FILE_MODE_WRITE );

      PetscViewerASCIIOpen( PETSC_COMM_WORLD , outname.c_str() , &viewer );
      VecView( U , viewer );
    
      PetscViewerDestroy( &viewer );
    
      app->dt_output_counter += 1;

    }

    // Look for new input parameters from driver
    if ( (event_list[i] == 0) && ( t < app->t_final ) ) {

      const std::string name     = "T_" + std::to_string( (app->dt_counter + 1) * app->dt_check ).substr(0,6) + ".dat";
      const std::string petscout = "Attepting to read new T(amplitude,x_mean,y_mean,sigma) at t = %5.4f seconds from file " + name + "\n";
      PetscPrintf( PETSC_COMM_WORLD , petscout.c_str() , (double)t );

      // Output to file that you are ready for a new value
      const std::string outname  = "complete_" + std::to_string( (app->dt_counter + 1) * app->dt_check ).substr(0,6) + ".dat";
      std::ofstream fout(outname);
      fout << "complete\n";
      fout.close();

      // Wait until receive a new value of m from an input file
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
        
          compute_new_temperature_profile( app , T_amp , T_x , T_y , T_sigma );
          
          app->dt_counter += 1;
        
          PetscPrintf( PETSC_COMM_WORLD , "Changing (T_amp, T_x, T_y, T_sigma) at t = %5.4f seconds to ( %5.4f , %5.4f , %5.4f , %5.4f)\n" , (double)t , (double)T_amp , (double)T_x , (double)T_y , (double)T_sigma );

          fin.close();        
          break;
        
        }

      }
      
      // Recompute ch parameters based on new temperature
      compute_eps2_and_sigma_from_temperature( ctx );
    }
    
  }

  return(0);

}

void compute_new_temperature_profile( AppCtx* user , PetscScalar T_amp , PetscScalar T_x , PetscScalar T_y , PetscScalar T_sigma  ) {

  DM             da   =user->da;
  PetscInt       i,j,xs,ys,xm,ym,Mx,My;
  PetscScalar    **T;
  PetscReal      hx,hy,x,y,r;
  
  PetscFunctionBeginUser;
  DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);

  hx = (user->Lx)/(PetscReal)(Mx-1);
  hy = (user->Ly)/(PetscReal)(My-1);

  /* Get pointers to vector data */
  DMDAVecGetArray(da,user->temperature,&T);

  /* Get local grid boundaries */
  DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);

  // Interior
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      T[j][i]     = T_amp * PetscExpReal( -0.5 * ( (j-T_x)*(j-T_x) + (i-T_y)*(i-T_y) ) / (T_sigma * T_sigma) );
    }
  }

  /* Restore vectors */
  DMDAVecRestoreArray(da,user->temperature,&T);
  
  return;


}
