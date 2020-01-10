#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <stdio.h>
#include "boundary_conditions.h"

PetscReal compute_residuals_no_explicit_boundary_resets( PetscReal*** uarray , PetscReal rhs_ijk , PetscReal udot_ijk , PetscInt Mx , PetscInt My , PetscInt Mz , PetscInt i , PetscInt j , PetscInt k ) {

  // Just compute residuals, no resets for boundaries required
  
  return udot_ijk - rhs_ijk;

}

void set_boundary_ghost_nodes_normal_extrapolation(  AppCtx* user , PetscScalar*** uarray , PetscInt Mx , PetscInt My , PetscInt Mz , PetscInt i , PetscInt j , PetscInt k ) {

  // Extrapolate out in the normal direction to fill ghost node layers
  // WARNING: this does NOT fill in any of the corner points and hence this method should
  // only be used for boundary conditions where a "star" stencil is needed for uarray (NOT a "box")
  
  if (i <= 1) {
    uarray[k][j][-1] = 2 * uarray[k][j][0] - uarray[k][j][1];
    uarray[k][j][-2] = 3 * uarray[k][j][0] - 2 * uarray[k][j][1];
  }
  
  else if (i >= Mx-2) {
    uarray[k][j][Mx]   = 2 * uarray[k][j][Mx-1] - uarray[k][j][Mx-2];
    uarray[k][j][Mx+1] = 3 * uarray[k][j][Mx-1] - 2 * uarray[k][j][Mx-2];
  }
  
  if (j <= 1) {
    uarray[k][-1][i] = 2 * uarray[k][0][i] - uarray[k][1][i];
    uarray[k][-2][i] = 3 * uarray[k][0][i] - 2 * uarray[k][1][i];    
  }
  
  else if (j >= My-2) {
    uarray[k][My][i]   = 2 * uarray[k][My-1][i] - uarray[k][My-2][i];
    uarray[k][My+1][i] = 3 * uarray[k][My-1][i] - 2 * uarray[k][My-2][i];
  }
  
  if (k <= 1) {
    uarray[-1][j][i] = 2 * uarray[0][j][i] - uarray[1][j][i];
    uarray[-2][j][i] = 3 * uarray[0][j][i] - 2 * uarray[1][j][i];
  }
  
  else if (k >= Mz-2) {
    uarray[Mz][j][i]   = 2 * uarray[Mz-1][j][i] - uarray[Mz-2][j][i];
    uarray[Mz+1][j][i] = 3 * uarray[Mz-1][j][i] - 2 * uarray[Mz-2][j][i];
  }

};

// void set_boundary_ghost_nodes( AppCtx* user , PetscScalar*** uarray , PetscInt Mx , PetscInt My , PetscInt Mz , PetscInt i , PetscInt j , PetscInt k ) {
  
//   if ( (user->boundary == 0) || (user->boundary == 3) || (user->boundary == 4) ) {
//     // 0: Dirichlet
//     // 3,4: mixed dirichlet-neumann, just fill ghost cells with dirichlet and we will reset boundary residuals later
//     if (i <= 1) {
//       for (int kk = k-2 ; kk <= k+2 ; kk++) {
//         for (int jj = j-2 ; jj <= j+2 ; jj++) {
//           uarray[kk][jj][-1] = user->dirichlet_bc;
//           uarray[kk][jj][-2] = user->dirichlet_bc;
//         }      
//       }
//     }
      
//     else if (i >= Mx-2) {
//       for (int kk = k-2 ; kk <= k+2 ; kk++) {
//         for (int jj = j-2 ; jj <= j+2 ; jj++) {
//           uarray[kk][jj][Mx]   = user->dirichlet_bc;
//           uarray[kk][jj][Mx+1] = user->dirichlet_bc;
//         }
//       }
//     }
    
//     if (j <= 1) {
//       for (int kk = k-2 ; kk <= k+2 ; kk++) {
//         for (int ii = i-2 ; ii <= i+2 ; ii++) {
//           uarray[kk][-2][ii] = user->dirichlet_bc;
//           uarray[kk][-1][ii] = user->dirichlet_bc;
//         }
//       }
//     }
    
//     else if (j >= My-2) {
//       for (int kk = k-2 ; kk <= k+2 ; kk++) {
//         for (int ii = i-2 ; ii <= i+2 ; ii++) {
//           uarray[kk][My][ii]   = user->dirichlet_bc;
//           uarray[kk][My+1][ii] = user->dirichlet_bc;
//         }
//       }
//     }
    
//     if (k <= 1) {
//       for (int jj = j-2 ; jj <= j+2 ; jj++) {
//         for (int ii = i-2 ; ii <= i+2 ; ii++) {
//           uarray[-1][jj][ii] = user->dirichlet_bc;
//           uarray[-2][jj][ii] = user->dirichlet_bc;
//         }      
//       }
//     }
    
//     else if (k >= Mz-2) {
//       for (int jj = j-2 ; jj <= j+2 ; jj++) {
//         for (int ii = i-2 ; ii <= i+2 ; ii++) {
//           uarray[Mz][jj][ii]   = user->dirichlet_bc;
//           uarray[Mz+1][jj][ii] = user->dirichlet_bc;
//         }
//       }
//     }

//   }

//   else if ( user->boundary == 1 ) {
//     // 1: Neumann
//     if (i <= 1) {
//       uarray[k][j][-2] = uarray[k][j][2];
//       uarray[k][j][-1] = uarray[k][j][1];
//     }
//     if (i >= Mx-2) {
//       uarray[k][j][Mx]   = uarray[k][j][Mx-2];
//       uarray[k][j][Mx+1] = uarray[k][j][Mx-3];
//     }
//     if (j <= 1) {
//       uarray[k][-2][i] = uarray[k][2][i];
//       uarray[k][-1][i] = uarray[k][1][i];
//     }
//     if (j >= My-2) {
//       uarray[k][My][i]   = uarray[k][My-2][i];
//       uarray[k][My+1][i] = uarray[k][My-3][i];
//     }
//     if (k <= 1) {
//       uarray[-2][j][i] = uarray[2][j][i];
//       uarray[-1][j][i] = uarray[1][j][i];
//     }
//     if (k >= Mz-2) {
//       uarray[Mz][j][i]   = uarray[Mz-2][j][i];
//       uarray[Mz+1][j][i] = uarray[Mz-3][j][i];
//     }
    
//   }

//   else if ( user->boundary == 2 ) {
//     // Periodic
//   }

//   return;
  
// };

PetscReal reset_boundary_residual_values_for_neumann_bc( PetscReal*** uarray , PetscReal rhs_ijk , PetscReal udot_ijk , PetscInt Mx , PetscInt My , PetscInt Mz , PetscInt i , PetscInt j , PetscInt k ) {

  PetscReal f_kji;
  
  if (i == 0 && j == 0 && k == 0) {               // SW-bottom corner
    f_kji = uarray[k][j][i] - uarray[k+1][j+1][i+1];
  }
  else if (i == 0 && j == 0 && k == Mz-1) {       // SW-top corner
    f_kji = uarray[k][j][i] - uarray[k-1][j+1][i+1];
  }
  else if (i == Mx-1 && j == 0 && k == 0) {          // SE-bottom corner 
    f_kji = uarray[k][j][i] - uarray[k+1][j+1][i-1];
  }
  else if (i == Mx-1 && j == 0 && k == Mz-1) {       // SE-top corner 
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
  
  else if ( i == 0 ) {        // W-face
    f_kji = uarray[k][j][i] - uarray[k][j][i+1];
  }
  else if ( i == Mx-1 ) {  // E-face
    f_kji = uarray[k][j][i] - uarray[k][j][i-1];
  }
  else if ( j == 0 ) {        // S-face 
    f_kji = uarray[k][j][i] - uarray[k][j+1][i];
  }
  else if ( j == My-1 ) {  // N-face
    f_kji = uarray[k][j][i] - uarray[k][j-1][i];
  }
  else if ( k == 0 ) {        // Bottom-face 
    f_kji = uarray[k][j][i] - uarray[k+1][j][i];
  }
  else if ( k == Mz-1 ) {  // Top-face
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

void set_boundary_ghost_nodes_dirichlet_singleframe( AppCtx* user , PetscScalar*** uarray , PetscInt Mx , PetscInt My , PetscInt Mz , PetscInt i , PetscInt j , PetscInt k ) {

  // Set three points around center in ghost regions

  if (i == 0) {
    for (int kk = k-1 ; kk <= k+1 ; kk++) {
      for (int jj = j-1 ; jj <= j+1 ; jj++) {
        uarray[kk][jj][-1] = user->dirichlet_bc;
      }
    }
  }

  else if (i == Mx-1) {
    for (int kk = k-1 ; kk <= k+1 ; kk++) {
      for (int jj = j-1 ; jj <= j+1 ; jj++) {
        uarray[kk][jj][Mx]   = user->dirichlet_bc;
      }
    }
  }
    
  if (j == 0) {
    for (int kk = k-1 ; kk <= k+1 ; kk++) {
      for (int ii = i-1 ; ii <= i+1 ; ii++) {
        uarray[kk][-1][ii] = user->dirichlet_bc;
      }
    }
  }
  
  else if (j == My-1) {
    for (int kk = k-1 ; kk <= k+1 ; kk++) {
      for (int ii = i-1 ; ii <= i+1 ; ii++) {
        uarray[kk][My][ii]   = user->dirichlet_bc;
      }
    }
  }  

  if (k == 0) {
    for (int jj = j-1 ; jj <= j+1 ; jj++) {
      for (int ii = i-1 ; ii <= i+1 ; ii++) {
        uarray[-1][jj][ii] = user->dirichlet_bc;
      }
    }
  }
  
  else if (k == Mz-1) {
    for (int jj = j-1 ; jj <= j+1 ; jj++) {
      for (int ii = i-1 ; ii <= i+1 ; ii++) {
        uarray[Mz][jj][ii]   = user->dirichlet_bc;
      }
    }
  }  
  
  return;
  
};
