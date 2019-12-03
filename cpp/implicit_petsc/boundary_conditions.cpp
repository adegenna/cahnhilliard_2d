#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <stdio.h>
#include "boundary_conditions.h"

void set_boundary_ghost_nodes_dirichlet( AppCtx* user , PetscScalar** uarray , PetscInt Mx , PetscInt My , PetscInt i , PetscInt j ) {

  // Set three points around center in ghost regions

  if (i <= 1) {
    for (int jj = j-2 ; jj <= j+2 ; jj++) {
      uarray[jj][-1] = user->dirichlet_bc;
      uarray[jj][-2] = user->dirichlet_bc;
    }
  }

  else if (i >= Mx-2) {
    for (int jj = j-2 ; jj <= j+2 ; jj++) {
      uarray[jj][Mx]   = user->dirichlet_bc;
      uarray[jj][Mx+1] = user->dirichlet_bc;
    }
  }
    
  if (j <= 1) {
    for (int ii = i-2 ; ii <= i+2 ; ii++) {
      uarray[-1][ii] = user->dirichlet_bc;
      uarray[-2][ii] = user->dirichlet_bc;
    }    
  }
  
  else if (j >= My-2) {
    for (int ii = i-2 ; ii <= i+2 ; ii++) {
      uarray[My][ii]   = user->dirichlet_bc;
      uarray[My+1][ii] = user->dirichlet_bc;
    }    
  }  
  
  return;
  
};

void set_boundary_ghost_nodes_neumann( AppCtx* user , PetscScalar** uarray , PetscInt Mx , PetscInt My , PetscInt i , PetscInt j ) {
  
  if (i <= 1) {
    uarray[j][-2] = uarray[j][2];
    uarray[j][-1] = uarray[j][1];
  }
  if (i >= Mx-2) {
    uarray[j][Mx]   = uarray[j][Mx-2];
    uarray[j][Mx+1] = uarray[j][Mx-3];
  }
  if (j <= 1) {
    uarray[-2][i] = uarray[2][i];
    uarray[-1][i] = uarray[1][i];
  }
  if (j >= My-2) {
    uarray[My][i]   = uarray[My-2][i];
    uarray[My+1][i] = uarray[My-3][i];
  }
  
  return;
  
};


ThirteenPointStencil get_thirteen_point_stencil( AppCtx* user , PetscReal** uarray , PetscInt Mx , PetscInt My , PetscInt i , PetscInt j ) {
  
  ThirteenPointStencil stencil;
  
  stencil.c_i   = uarray[j][i];
  stencil.c_jm2 = uarray[j-2][i];
  stencil.c_jm1 = uarray[j-1][i];
  stencil.c_jp1 = uarray[j+1][i];
  stencil.c_jp2 = uarray[j+2][i];
  stencil.c_im2 = uarray[j][i-2];
  stencil.c_im1 = uarray[j][i-1];
  stencil.c_ip1 = uarray[j][i+1];
  stencil.c_ip2 = uarray[j][i+2];
  stencil.c_ul  = uarray[j-1][i-1];
  stencil.c_ur  = uarray[j-1][i+1];
  stencil.c_bl  = uarray[j+1][i-1];
  stencil.c_br  = uarray[j+1][i+1];

  return stencil;
  
}

PetscReal reset_boundary_residual_values_for_neumann_bc( PetscReal** uarray , PetscReal rhs_ij , PetscReal udot_ij , PetscInt Mx , PetscInt My , PetscInt i , PetscInt j ) {

  PetscReal f_ji;
  
  if (i == 0 && j == 0) {            // SW corner
    f_ji = uarray[j][i] - uarray[j+1][i+1];
  }
  else if (i == Mx-1 && j == 0) {    // SE corner 
    f_ji = uarray[j][i] - uarray[j+1][i-1];
  }
  else if (i == 0 && j == My-1) {    // NW corner 
    f_ji = uarray[j][i] - uarray[j-1][i+1];
  }
  else if (i == Mx-1 && j == My-1) { // NE corner 
    f_ji = uarray[j][i] - uarray[j-1][i-1];
  }
	  
  else if ( (i == 0) || (i == 1) ) {        // Left 
    f_ji = uarray[j][i] - uarray[j][i+1];
  }
  else if ( (i == Mx-1) || (i == Mx-2) ) {  // Right 
    f_ji = uarray[j][i] - uarray[j][i-1];
  }
  else if ( (j == 0) || (j == 1) ) {        // Bottom 
    f_ji = uarray[j][i] - uarray[j+1][i];
  }
  else if ( (j == My-1) || (j == My-2) ) {  // Top 
    f_ji = uarray[j][i] - uarray[j-1][i];
  }
  else
    f_ji = udot_ij - rhs_ij;


  return f_ji;
}

PetscReal reset_boundary_residual_values_for_dirichlet_bottom_neumann_remainder_bc( PetscReal** uarray , PetscReal rhs_ij , PetscReal udot_ij , PetscInt Mx , PetscInt My , PetscInt i , PetscInt j ) {

  PetscReal f_ji;

  // Default: interior value
  f_ji = udot_ij - rhs_ij;

  // Neumann parts
  if (i == Mx-1 && j == 0) {    // SE corner 
    f_ji = uarray[j][i] - uarray[j+1][i-1];
  }
  else if (i == Mx-1 && j == My-1) { // NE corner 
    f_ji = uarray[j][i] - uarray[j-1][i-1];
  }	  
  else if ( (i == Mx-1) || (i == Mx-2) ) {  // Right 
    f_ji = uarray[j][i] - uarray[j][i-1];
  }
  else if ( (j == 0) || (j == 1) ) {        // Bottom 
    f_ji = uarray[j][i] - uarray[j+1][i];
  }
  else if ( (j == My-1) || (j == My-2) ) {  // Top 
    f_ji = uarray[j][i] - uarray[j-1][i];
  }

  // Dirichlet parts
  if ( (i == 0) || (i == 1) ) {        // Left 
    f_ji = udot_ij - rhs_ij;
  }

  return f_ji;
}

PetscReal reset_boundary_residual_values_for_dirichlet_topandbottom_neumann_remainder_bc( PetscReal** uarray , PetscReal rhs_ij , PetscReal udot_ij , PetscInt Mx , PetscInt My , PetscInt i , PetscInt j ) {

  PetscReal f_ji;

  // Default: interior value
  f_ji = udot_ij - rhs_ij;

  // Neumann parts
  if ( (j == 0) || (j == 1) ) {        // Bottom 
    f_ji = uarray[j][i] - uarray[j+1][i];
  }
  else if ( (j == My-1) || (j == My-2) ) {  // Top 
    f_ji = uarray[j][i] - uarray[j-1][i];
  }

  // Dirichlet parts
  if ( (i == 0) || (i == 1) ) {        // Left 
    f_ji = udot_ij - rhs_ij;
  }
  else if ( (i == Mx-1) || (i == Mx-2) ) { // Right
    f_ji = udot_ij - rhs_ij;
  }

  return f_ji;
}

PetscReal compute_residuals_no_explicit_boundary_resets( PetscReal** uarray , PetscReal rhs_ij , PetscReal udot_ij , PetscInt Mx , PetscInt My , PetscInt i , PetscInt j ) {

  // Just compute residuals, no resets for boundaries required
  
  return udot_ij - rhs_ij;

}



ThirteenPointStencil apply_dirichlet_bc( AppCtx* user , PetscReal** uarray , PetscInt Mx , PetscInt My , PetscInt i , PetscInt j ) {

  ThirteenPointStencil stencil;

  // Computes 13-pt stencil for dirichlet bcs

  if ((j-2 >= 0) && (j-2 <= My-1) && (i >= 0) && (i <= Mx-1)  )    //c_jm2 
    stencil.c_jm2 = uarray[j-2][i];
  else
    stencil.c_jm2 = user->dirichlet_bc;

  if ((j-1 >= 0) && (j-1 <= My-1) && (i >= 0) && (i <= Mx-1)  )    //c_jm1
    stencil.c_jm1 = uarray[j-1][i];
  else
    stencil.c_jm1 = user->dirichlet_bc;

  if ((j+1 >= 0) && (j+1 <= My-1) && (i >= 0) && (i <= Mx-1)  )    //c_jp1
    stencil.c_jp1 = uarray[j+1][i];
  else
    stencil.c_jp1 = user->dirichlet_bc;

  if ((j+2 >= 0) && (j+2 <= My-1) && (i >= 0) && (i <= Mx-1)  )    //c_jp2
    stencil.c_jp2 = uarray[j+2][i];
  else
    stencil.c_jp2 = user->dirichlet_bc;

  if ((j >= 0) && (j <= My-1) && (i-2 >= 0) && (i-2 <= Mx-1)  )    //c_im2 
    stencil.c_im2 = uarray[j][i-2];
  else
    stencil.c_im2 = user->dirichlet_bc;

  if ((j >= 0) && (j <= My-1) && (i-1 >= 0) && (i-1 <= Mx-1)  )    //c_im1 
    stencil.c_im1 = uarray[j][i-1];
  else
    stencil.c_im1 = user->dirichlet_bc;

  if ((j >= 0) && (j <= My-1) && (i+1 >= 0) && (i+1 <= Mx-1)  )    //c_ip1 
    stencil.c_ip1 = uarray[j][i+1];
  else
    stencil.c_ip1 = user->dirichlet_bc;

  if ((j >= 0) && (j <= My-1) && (i+2 >= 0) && (i+2 <= Mx-1)  )    //c_ip2 
    stencil.c_ip2 = uarray[j][i+2];
  else
    stencil.c_ip2 = user->dirichlet_bc;

  if ((j-1 >= 0) && (j-1 <= My-1) && (i-1 >= 0) && (i-1 <= Mx-1)  )    //c_ul
    stencil.c_ul = uarray[j-1][i-1];
  else
    stencil.c_ul = user->dirichlet_bc;

  if ((j-1 >= 0) && (j-1 <= My-1) && (i+1 >= 0) && (i+1 <= Mx-1)  )    //c_ur
    stencil.c_ur = uarray[j-1][i+1];
  else
    stencil.c_ur = user->dirichlet_bc;

  if ((j+1 >= 0) && (j+1 <= My-1) && (i-1 >= 0) && (i-1 <= Mx-1)  )    //c_bl
    stencil.c_bl = uarray[j+1][i-1];
  else
    stencil.c_bl = user->dirichlet_bc;

  if ((j+1 >= 0) && (j+1 <= My-1) && (i+1 >= 0) && (i+1 <= Mx-1)  )    //c_br
    stencil.c_br = uarray[j+1][i+1];
  else
    stencil.c_br = user->dirichlet_bc;
  
  return stencil;

};



ThirteenPointStencil apply_neumann_bc( AppCtx* user , PetscReal** uarray , PetscInt Mx , PetscInt My , PetscInt i , PetscInt j ) {

  ThirteenPointStencil stencil;

  // Computes 13-pt stencil for neumann bcs

  stencil.c_i   = uarray[j][i];
  
  if ((j-2 >= 0) && (j-2 <= My-1) && (i >= 0) && (i <= Mx-1)  )    //c_jm2 
    stencil.c_jm2 = uarray[j-2][i];
  else
    stencil.c_jm2 = uarray[j+2][i];

  if ((j-1 >= 0) && (j-1 <= My-1) && (i >= 0) && (i <= Mx-1)  )    //c_jm1
    stencil.c_jm1 = uarray[j-1][i];
  else
    stencil.c_jm1 = uarray[j+1][i];

  if ((j+1 >= 0) && (j+1 <= My-1) && (i >= 0) && (i <= Mx-1)  )    //c_jp1
    stencil.c_jp1 = uarray[j+1][i];
  else
    stencil.c_jp1 = uarray[j-1][i];

  if ((j+2 >= 0) && (j+2 <= My-1) && (i >= 0) && (i <= Mx-1)  )    //c_jp2
    stencil.c_jp2 = uarray[j+2][i];
  else
    stencil.c_jp2 = uarray[j-2][i];

  if ((j >= 0) && (j <= My-1) && (i-2 >= 0) && (i-2 <= Mx-1)  )    //c_im2 
    stencil.c_im2 = uarray[j][i-2];
  else
    stencil.c_im2 = uarray[j][i+2];

  if ((j >= 0) && (j <= My-1) && (i-1 >= 0) && (i-1 <= Mx-1)  )    //c_im1 
    stencil.c_im1 = uarray[j][i-1];
  else
    stencil.c_im1 = uarray[j][i+1];

  if ((j >= 0) && (j <= My-1) && (i+1 >= 0) && (i+1 <= Mx-1)  )    //c_ip1 
    stencil.c_ip1 = uarray[j][i+1];
  else
    stencil.c_ip1 = uarray[j][i-1];

  if ((j >= 0) && (j <= My-1) && (i+2 >= 0) && (i+2 <= Mx-1)  )    //c_ip2 
    stencil.c_ip2 = uarray[j][i+2];
  else
    stencil.c_ip2 = uarray[j][i-2];

  if ( (j-1 >= 0) && (j-1 <= My-1) && (i-1 >= 0) && (i-1 <= Mx-1) )    //c_ul
    stencil.c_ul = uarray[j-1][i-1];
  else {
    if ( (j-1 >= 0) && (j-1 <= My-1) && (i+1 >= 0) && (i+1 <= Mx-1) )
      stencil.c_ul = uarray[j-1][i+1];
    else if ( (j+1 >= 0) && (j+1 <= My-1) && (i-1 >= 0) && (i-1 <= Mx-1) )
      stencil.c_ul = uarray[j+1][i-1];
    else
      stencil.c_ul = uarray[j+1][i+1];
  }
  
  if ( (j-1 >= 0) && (j-1 <= My-1) && (i+1 >= 0) && (i+1 <= Mx-1) )    //c_ur
    stencil.c_ur = uarray[j-1][i+1];
  else {
    if ( (j-1 >= 0) && (j-1 <= My-1) && (i-1 >= 0) && (i-1 <= Mx-1) )
      stencil.c_ur = uarray[j-1][i-1];
    else if ( (j+1 >= 0) && (j+1 <= My-1) && (i+1 >= 0) && (i+1 <= Mx-1) )
      stencil.c_ur = uarray[j+1][i+1];
    else
      stencil.c_ur = uarray[j+1][i-1];
  }
  
  if ((j+1 >= 0) && (j+1 <= My-1) && (i-1 >= 0) && (i-1 <= Mx-1)  )    //c_bl
    stencil.c_bl = uarray[j+1][i-1];
  else {
    if ( (j+1 >= 0) && (j+1 <= My-1) && (i+1 >= 0) && (i+1 <= Mx-1) )
      stencil.c_bl = uarray[j+1][i+1];
    else if ( (j-1 >= 0)  && (j-1 <= My-1) && (i-1 >= 0) && (i-1 <= Mx-1) )
      stencil.c_bl = uarray[j-1][i-1];
    else
      stencil.c_bl = uarray[j-1][i+1];
  }
  
  if ((j+1 >= 0) && (j+1 <= My-1) && (i+1 >= 0) && (i+1 <= Mx-1)  )    //c_br
    stencil.c_br = uarray[j+1][i+1];
  else {
    if ( (j+1 >= 0) && (j+1 <= My-1) && (i-1 >= 0) && (i-1 <= Mx-1) )
      stencil.c_br = uarray[j+1][i-1];
    else if ( (j-1 >= 0)  && (j-1 <= My-1) && (i+1 >= 0) && (i+1 <= Mx-1) )
      stencil.c_br = uarray[j-1][i+1];
    else
      stencil.c_br = uarray[j-1][i-1];
  }
  
  return stencil;

};
