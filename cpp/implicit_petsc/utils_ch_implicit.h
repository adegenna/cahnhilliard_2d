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
  DM        da;
  PetscReal c;
  PetscInt  boundary;            /* Type of boundary condition */
  PetscBool viewJacobian;
  PetscReal Lx, Ly;            // Length of domain in each direction
  PetscReal t_final;           // Final time of simulation
  PetscReal dirichlet_bc;      // Value of dirichlet bc
  PetscReal dt_check   = 0.1;  // Value of time increment where you change the parameters/temperature
  PetscInt  dt_counter = 0;    // Counter that keeps track of how many dt_check have gone by so far
  PetscReal dt_output  = 0.02; // Value of time increment where you change the parameters/temperature
  PetscInt  dt_output_counter = 0;   // Counter that keeps track of how many dt_output have gone by so far
  
  PetscReal m          = -0.15;  // CH parameter: value of m (avg concentration)
} AppCtx;

AppCtx parse_petsc_options( );

#endif
