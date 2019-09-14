#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <stdio.h>
#include "boundary_conditions.h"

void set_boundary_ghost_nodes( AppCtx* user , PetscScalar** uarray , PetscInt Mx , PetscInt My , PetscInt Mz , PetscInt i , PetscInt j , PetscInt k ) {
  
  if ( (user->boundary == 0) || (user->boundary == 3) || (user->boundary == 4) ) {
    // 0: Dirichlet
    // 3,4: mixed dirichlet-neumann, just fill ghost cells with dirichlet and we will reset boundary residuals later
    if (i <= 1) {
      uarray[k][j][-2] = user->dirichlet_bc;
      uarray[k][j][-1] = user->dirichlet_bc;
    }
    if (i >= Mx-2) {
      uarray[k][j][Mx]   = user->dirichlet_bc;
      uarray[k][j][Mx+1] = user->dirichlet_bc;
    }
    if (j <= 1) {
      uarray[k][-2][i] = user->dirichlet_bc;
      uarray[k][-1][i] = user->dirichlet_bc;
    }
    if (j >= My-2) {
      uarray[k][My][i]   = user->dirichlet_bc;
      uarray[k][My+1][i] = user->dirichlet_bc;
    }
    if (k <= 1) {
      uarray[-1][j][i] = user->dirichlet_bc;
      uarray[-2][j][i] = user->dirichlet_bc;
    }
    if (k >= Mz-2) {
      uarray[Mz][j][i]   = user->dirichlet_bc;
      uarray[Mz+1][j][i] = user->dirichlet_bc;
    }

  }

  else if ( user->boundary == 1 ) {
    // 1: Neumann
    if (i <= 1) {
      uarray[k][j][-2] = uarray[k][j][2];
      uarray[k][j][-1] = uarray[k][j][1];
    }
    if (i >= Mx-2) {
      uarray[k][j][Mx]   = uarray[k][j][Mx-2];
      uarray[k][j][Mx+1] = uarray[k][j][Mx-3];
    }
    if (j <= 1) {
      uarray[k][-2][i] = uarray[k][2][i];
      uarray[k][-1][i] = uarray[k][1][i];
    }
    if (j >= My-2) {
      uarray[k][My][i]   = uarray[k][My-2][i];
      uarray[k][My+1][i] = uarray[k][My-3][i];
    }
    if (k <= 1) {
      uarray[-2][j][i] = uarray[2][j][i];
      uarray[-1][j][i] = uarray[1][j][i];
    }
    if (k >= Mz-2) {
      uarray[Mz][j][i]   = uarray[Mz-2][j][i];
      uarray[Mz+1][j][i] = uarray[Mz-3][j][i];
    }
    
  }

  else if ( user->boundary == 2 ) {
    // Periodic
  }

  return;
  
};

PetscReal reset_boundary_residual_values_for_neumann_bc( PetscReal** uarray , PetscReal rhs_ijk , PetscReal udot_ijk , PetscInt Mx , PetscInt My , PetscInt Mz , PetscInt i , PetscInt j , PetscInt k ) {

  PetscReal f_kji;
  
  if (i == 0 && j == 0 && k == 0) {               // SW-bottom corner
    f_kji = uarray[k][j][i] - uarray[k+1][j+1][i+1];
  }
  else if (i == 0 && j == 0 && k == Mz-1) {       // SW-top corner
    f_kji = uarray[k][j][i] - uarray[k-1][j+1][i+1];
  }
  else if (i == Mx-1 && j == 0 && k=0) {          // SE-bottom corner 
    f_kji = uarray[k][j][i] - uarray[k+1][j+1][i-1];
  }
  else if (i == Mx-1 && j == 0 && k=Mz-1) {       // SE-top corner 
    f_kji = uarray[k][j][i] - uarray[k-1][j+1][i-1];
  }
  else if (i == 0 && j == My-1 && k == 0) {       // NW-bottom corner 
    f_kji = uarray[k][j][i] - uarray[k+1][j-1][i+1];
  }
  else if (i == 0 && j == My-1 && k == Mz-1) {    // NW-top corner 
    f_kji = uarray[k][j][i] - uarray[k-1][j-1][i+1];
  }
  else if (i == Mx-1 && j == My-1 && k == 0) {    // NE-bottom corner 
    f_kji = uarray[k][j][i] - uarray[k+1][j-1][i-1];
  }
  else if (i == Mx-1 && j == My-1 && k == Mz-1) { // NE-top corner 
    f_kji = uarray[k][j][i] - uarray[k-1][j-1][i-1];
  }
  
  else if ( (i == 0) || (i == 1) ) {        // W-face
    f_kji = uarray[k][j][i] - uarray[k][j][i+1];
  }
  else if ( (i == Mx-1) || (i == Mx-2) ) {  // E-face
    f_kji = uarray[k][j][i] - uarray[k][j][i-1];
  }
  else if ( (j == 0) || (j == 1) ) {        // S-face 
    f_kji = uarray[k][j][i] - uarray[k][j+1][i];
  }
  else if ( (j == My-1) || (j == My-2) ) {  // N-face
    f_kji = uarray[k][j][i] - uarray[k][j-1][i];
  }
  else if ( (k == 0) || (k == 1) ) {        // Bottom-face 
    f_kji = uarray[k][j][i] - uarray[k+1][j][i];
  }
  else if ( (k == Mz-1) || (k == Mz-2) ) {  // Top-face
    f_kji = uarray[k][j][i] - uarray[k-1][j][i];
  }
  
  else
    f_kji = udot_ijk - rhs_ijk;


  return f_kji;
}

PetscReal reset_boundary_residual_values_for_dirichlet_bottom_neumann_remainder_bc( PetscReal** uarray , PetscReal rhs_ijk , PetscReal udot_ij , PetscInt Mx , PetscInt My , PetscInt i , PetscInt j ) {

  PetscReal f_kji;

  // Default: interior value
  f_kji = udot_ij - rhs_ijk;

  // Neumann parts
  if (i == Mx-1 && j == 0) {    // SE corner 
    f_kji = uarray[j][i] - uarray[j+1][i-1];
  }
  else if (i == Mx-1 && j == My-1) { // NE corner 
    f_kji = uarray[j][i] - uarray[j-1][i-1];
  }	  
  else if ( (i == Mx-1) || (i == Mx-2) ) {  // Right 
    f_kji = uarray[j][i] - uarray[j][i-1];
  }
  else if ( (j == 0) || (j == 1) ) {        // Bottom 
    f_kji = uarray[j][i] - uarray[j+1][i];
  }
  else if ( (j == My-1) || (j == My-2) ) {  // Top 
    f_kji = uarray[j][i] - uarray[j-1][i];
  }

  // Dirichlet parts
  if ( (i == 0) || (i == 1) ) {        // Left 
    f_kji = udot_ij - rhs_ijk;
  }

  return f_kji;
}

PetscReal reset_boundary_residual_values_for_dirichlet_topandbottom_neumann_remainder_bc( PetscReal** uarray , PetscReal rhs_ijk , PetscReal udot_ij , PetscInt Mx , PetscInt My , PetscInt i , PetscInt j ) {

  PetscReal f_kji;

  // Default: interior value
  f_kji = udot_ij - rhs_ijk;

  // Neumann parts
  if ( (j == 0) || (j == 1) ) {        // Bottom 
    f_kji = uarray[j][i] - uarray[j+1][i];
  }
  else if ( (j == My-1) || (j == My-2) ) {  // Top 
    f_kji = uarray[j][i] - uarray[j-1][i];
  }

  // Dirichlet parts
  if ( (i == 0) || (i == 1) ) {        // Left 
    f_kji = udot_ij - rhs_ijk;
  }
  else if ( (i == Mx-1) || (i == Mx-2) ) { // Right
    f_kji = udot_ij - rhs_ijk;
  }

  return f_kji;
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
