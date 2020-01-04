#ifndef __UTILS_CH_IMPLICIT_H__
#define __UTILS_CH_IMPLICIT_H__

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <stdio.h>

class ThirteenPointStencil {

  public:

    ThirteenPointStencil() { };
   ~ThirteenPointStencil() { };
    PetscReal c_i,c_im2,c_im1,c_ip1,c_ip2,c_jm2,c_jm1,c_jp1,c_jp2, c_ul,c_ur,c_bl,c_br;

};

/* AppCtx: used by FormIFunction() and FormIJacobian() */
typedef struct {
  DM        da_c , da_phi , da_T , pack;
  std::string physics = "ch";       // "ch": CH-only ; "thermal": thermal diffusion only ; "coupled_ch_thermal": coupled thermal-CH solver
  PetscReal c;

  std::string boundary  = "neumann";
  PetscReal (*residualFunction)( PetscReal** , PetscReal , PetscReal , PetscInt , PetscInt , PetscInt , PetscInt );
  
  PetscReal Lx, Ly;            // Length of domain in each direction
  PetscReal t_final;           // Final time of simulation
  PetscReal dirichlet_bc;      // Value of dirichlet bc
  PetscReal dt_check;          // Value of time increment where you change the parameters/temperature
  PetscInt  dt_counter = 0;    // Counter that keeps track of how many dt_check have gone by so far
  PetscReal dt_output;         // Value of time increment where you change the parameters/temperature
  PetscInt  dt_output_counter = 0;   // Counter that keeps track of how many dt_output have gone by so far
  std::string time_stepper    = "implicit";   // "implicit" or "explicit"
  PetscScalar dt          = 0.005;   // Default dt

  // Polymer physics defaults
  PetscScalar X_min    = 0.055;
  PetscScalar X_max    = 0.5;
  PetscScalar N        = 0.5 * ( 200.0 + 2000.0 );
  PetscScalar L_repeat = 0.5 * ( 20.0 + 80.0 );   // nanometers
  PetscScalar n_repeat = 15.0;
  PetscScalar L_omega  = n_repeat * L_repeat;
  PetscScalar L_kuhn   = 0.5 * ( 0.5 + 3.0 );     // nanometers

  // Thermal dynamics defaults
  PetscScalar T_min     = 0.1;
  PetscScalar T_max     = 1.0;
  std::string initial_temperature_file        = "initial_temperature.dat"; // File that holds the initial temperature field
  std::string initial_temperature_source_file = "initial_temperature_source.dat"; // File that holds the initial temperature source field
  std::string initial_soln_file               = "initial_soln.dat"; // File that holds the initial solution field
  PetscScalar D_T       = 1.0;  // Thermal diffusion coefficient
  
  // CH paramater defaults
  PetscScalar m         = 0.1;  // CH parameter: value of m (avg concentration)
  PetscScalar eps2_min  = 0.0;
  PetscScalar eps2_max  = 1.0;
  PetscScalar sigma_min = 0.0;
  PetscScalar sigma_max = pow( 10.0 , 10 );
  Vec         eps_2,sigma,temperature_source,X;
  
} AppCtx;

AppCtx parse_petsc_options( );

#endif
