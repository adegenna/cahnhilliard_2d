#include "utils_ch.h"

void compute_ch_nonlocal(const std::vector<double> &c,
			 std::vector<double> &dcdt,
			 const double t,
			 CHparamsVector& chpV,
			 SimInfo& info) {

  // Computes deterministic nonlocal CH dynamics
  // dc/dt = laplacian( u*c^3 - b*c ) - eps_2*biharm(c) - sigma*(c - m)

  // evaluate the second order term, 5 point central stencil
  # pragma omp parallel for
  for (int i = 0; i < info.ny; ++i) {
    for (int j = 0; j < info.nx; ++j) {
        
      const double c_i   = laplace_component( info.idx2d(i, j)      , c , chpV.u , chpV.b );
      const double c_im1 = laplace_component( info.idx2d(i - 1, j)  , c , chpV.u , chpV.b );
      const double c_ip1 = laplace_component( info.idx2d(i + 1, j)  , c , chpV.u , chpV.b );
      const double c_jm1 = laplace_component( info.idx2d(i, j - 1)  , c , chpV.u , chpV.b );
      const double c_jp1 = laplace_component( info.idx2d(i, j + 1)  , c , chpV.u , chpV.b );

      double dxx = (1.0 / (info.dx * info.dx)) * ( c_jp1 + c_jm1 - 2.0 * c_i );
      double dyy = (1.0 / (info.dy * info.dy)) * ( c_ip1 + c_im1 - 2.0 * c_i );
      dcdt[info.idx2d(i, j)] = dxx + dyy;
      
    }
  }

  // evaluate the 4th order term, 9 point central stencil
  # pragma omp parallel for
  for (int i = 0; i < info.ny; ++i) {
    for (int j = 0; j < info.nx; ++j) {
      
      const double c_i   = c[info.idx2d(i, j)];
      const double c_im1 = c[info.idx2d(i - 1, j)];
      const double c_ip1 = c[info.idx2d(i + 1, j)];
      const double c_im2 = c[info.idx2d(i - 2, j)];
      const double c_ip2 = c[info.idx2d(i + 2, j)];
      const double c_jm1 = c[info.idx2d(i, j - 1)];
      const double c_jp1 = c[info.idx2d(i, j + 1)];
      const double c_jm2 = c[info.idx2d(i, j - 2)];
      const double c_jp2 = c[info.idx2d(i, j + 2)];
      const double c_ul  = c[info.idx2d(i-1 , j-1)];
      const double c_ur  = c[info.idx2d(i-1 , j+1)];
      const double c_bl  = c[info.idx2d(i+1 , j-1)];
      const double c_br  = c[info.idx2d(i+1 , j+1)];

      // y-direction u_yyyy
      double dyyyy = 1.0 / (info.dy * info.dy * info.dy * info.dy) * 
        (c_ip2 - 4.0*c_ip1 + 6.0*c_i - 4.0*c_im1 + c_im2);

      // x-direction u_xxxx
      double dxxxx = 1.0 / (info.dx * info.dx * info.dx * info.dx) * 
        (c_jp2 - 4.0*c_jp1 + 6.0*c_i - 4.0*c_jm1 + c_jm2);

      // mixed term 2*u_xxyy
      double dxxyy = 1.0 / (info.dx * info.dx * info.dy * info.dy) * 
        2 * (4*c_i - 2*(c_im1 + c_ip1 + c_jm1 + c_jp1) + c_ul + c_ur + c_bl + c_br );

      dcdt[info.idx2d(i,j)] += -chpV.eps_2[info.idx2d(i,j)] * ( dxxxx + dyyyy + dxxyy );
      
    }
  }

  // evaluate linear term
  # pragma omp parallel for
  for (int i = 0; i < info.ny; ++i){
    for (int j = 0; j < info.nx; ++j){
        
      const double c_i        = c[info.idx2d(i, j)];
      dcdt[info.idx2d(i,j)]  += -chpV.sigma[info.idx2d(i,j)] * ( c_i - chpV.m[info.idx2d(i,j)] );

    }
  }


}

std::vector<double>& apply_dirichlet_bc( std::vector<double>& c ,
					 SimInfo& info ) {

  # pragma omp parallel for
  for (int i = 0; i < info.nx; ++i) {
    
    c[info.idx2d(0, i)]         = info.BC_dirichlet_ch;
    c[info.idx2d(1, i)]         = info.BC_dirichlet_ch;
    c[info.idx2d(info.ny-1, i)] = info.BC_dirichlet_ch;
    c[info.idx2d(info.ny-2, i)] = info.BC_dirichlet_ch;

  }

  # pragma omp parallel for
  for (int i = 0; i < info.ny; ++i) {

    c[info.idx2d(i, 0)]         = info.BC_dirichlet_ch;
    c[info.idx2d(i, 1)]         = info.BC_dirichlet_ch;
    c[info.idx2d(i, info.nx-1)] = info.BC_dirichlet_ch;
    c[info.idx2d(i, info.nx-2)] = info.BC_dirichlet_ch;
    
  }

  return c;

}

std::vector<double>& apply_neumann_bc( std::vector<double>& c ,
				       SimInfo& info ) {

  # pragma omp parallel for
  for (int i = 0; i < info.nx; ++i) {
    
    c[info.idx2d(0, i)]         = c[info.idx2d(4, i)];
    c[info.idx2d(1, i)]         = c[info.idx2d(3, i)];
    c[info.idx2d(info.ny-1, i)] = c[info.idx2d(info.ny-5, i)];
    c[info.idx2d(info.ny-2, i)] = c[info.idx2d(info.ny-4, i)];

  }
  
  # pragma omp parallel for
  for (int i = 0; i < info.ny; ++i) {

    c[info.idx2d(i, 0)]         = c[info.idx2d(i, 4)];
    c[info.idx2d(i, 1)]         = c[info.idx2d(i, 3)];
    c[info.idx2d(i, info.nx-1)] = c[info.idx2d(i, info.nx-5)];
    c[info.idx2d(i, info.nx-2)] = c[info.idx2d(i, info.nx-4)];

  }
  
  return c;

}

std::vector<double>& apply_mixed_bc_neumann_with_top_dirichlet( std::vector<double>& c ,
                                                                SimInfo& info ) {

  c = apply_neumann_bc( c , info );
  
  # pragma omp parallel for
  for (int i = 0; i < info.nx; ++i) {

    c[info.idx2d(info.ny-2, i)] = info.BC_dirichlet_ch;
    c[info.idx2d(info.ny-1, i)] = info.BC_dirichlet_ch;

  }

  return c;
  
}

std::vector<double>& apply_mixed_bc_neumann_with_bottom_dirichlet( std::vector<double>& c ,
								   SimInfo& info ) {

  c = apply_neumann_bc( c , info );
  
  # pragma omp parallel for
  for (int i = 0; i < info.nx; ++i) {
    
    c[info.idx2d(0, i)] = info.BC_dirichlet_ch;
    c[info.idx2d(1, i)] = info.BC_dirichlet_ch;

  }

  return c;
  
}

std::vector<double>& set_boundary_values_to_zero( std::vector<double> &dcdt ,
						  SimInfo& info ) {

  # pragma omp parallel for
  for (int i = 0; i < info.nx; ++i) {
    
    dcdt[info.idx2d(0, i)]         = 0;
    dcdt[info.idx2d(1, i)]         = 0;
    dcdt[info.idx2d(info.ny-1, i)] = 0;
    dcdt[info.idx2d(info.ny-2, i)] = 0;

  }

  # pragma omp parallel for
  for (int i = 0; i < info.ny; ++i) {

    dcdt[info.idx2d(i, 0)]         = 0;
    dcdt[info.idx2d(i, 1)]         = 0;
    dcdt[info.idx2d(i, info.nx-1)] = 0;
    dcdt[info.idx2d(i, info.nx-2)] = 0;

  }

  return dcdt;
  
}

void compute_ch_nonlocal_stationary_boundaries(const std::vector<double> &c,
					       std::vector<double> &dcdt,
					       const double t,
					       CHparamsVector& chpV,
					       SimInfo& info) {

  compute_ch_nonlocal(c, dcdt, t, chpV, info);
  dcdt = set_boundary_values_to_zero( dcdt , info );

}

std::vector<double>& freeze_corners( std::vector<double>& dcdt , SimInfo& info ) {

  dcdt[info.idx2d(0,0)] = 0;         dcdt[info.idx2d(0,1)] = 0;         dcdt[info.idx2d(0,info.nx-2)] = 0;         dcdt[info.idx2d(0,info.nx-1)] = 0;
  dcdt[info.idx2d(1,0)] = 0;         dcdt[info.idx2d(1,1)] = 0;         dcdt[info.idx2d(1,info.nx-2)] = 0;         dcdt[info.idx2d(1,info.nx-1)] = 0;
  dcdt[info.idx2d(info.ny-2,0)] = 0; dcdt[info.idx2d(info.ny-2,1)] = 0; dcdt[info.idx2d(info.ny-2,info.nx-2)] = 0; dcdt[info.idx2d(info.ny-2,info.nx-1)] = 0;
  dcdt[info.idx2d(info.ny-1,0)] = 0; dcdt[info.idx2d(info.ny-1,1)] = 0; dcdt[info.idx2d(info.ny-1,info.nx-2)] = 0; dcdt[info.idx2d(info.ny-1,info.nx-1)] = 0;

  return dcdt;
  
}

void compute_ch_nonlocal_neumannBC(const std::vector<double> &c,
                                   std::vector<double> &dcdt,
                                   const double t,
                                   CHparamsVector& chpV,
                                   SimInfo& info) {

  compute_ch_nonlocal(c, dcdt, t, chpV, info);
  dcdt = apply_neumann_bc( dcdt , info );
  dcdt = freeze_corners( dcdt , info ); // 4 corners don't affect dynamics
  
}

void compute_ch_nonlocal_mixedBC_neumann_with_bottom_dirichlet(const std::vector<double> &c,
                                                               std::vector<double> &dcdt,
                                                               const double t,
                                                               CHparamsVector& chpV,
                                                               SimInfo& info) {

  compute_ch_nonlocal(c, dcdt, t, chpV, info);
  dcdt = apply_mixed_bc_neumann_with_bottom_dirichlet( dcdt , info );
  dcdt = freeze_corners( dcdt , info );
  
}

void compute_ch_nonlocal_mixedBC_neumann_with_top_dirichlet(const std::vector<double> &c,
                                                            std::vector<double> &dcdt,
                                                            const double t,
                                                            CHparamsVector& chpV,
                                                            SimInfo& info) {

  compute_ch_nonlocal(c, dcdt, t, chpV, info);
  dcdt = apply_mixed_bc_neumann_with_top_dirichlet( dcdt , info );
  dcdt = freeze_corners( dcdt , info );
  
}

double laplace_component(int i ,
                         const std::vector<double>& c ,
                         const std::vector<double>& u ,
                         const std::vector<double>& b ) {

  return u[i] * (c[i] * c[i] * c[i]) - b[i] * c[i];
}

CHparamsVector compute_chparams_using_temperature( CHparamsVector& chpV0,
                                                   SimInfo& info,
                                                   std::vector<double> T ) {

  CHparamsVector chpV = chpV0;
  double deps2_dT     = ( chpV.eps2_max  - chpV.eps2_min )  / ( chpV.T_max - chpV.T_min );
  double dsigma_dT    = ( chpV.sigma_max - chpV.sigma_min ) / ( chpV.T_max - chpV.T_min );
  
  # pragma omp parallel for
  for (int i = 0; i < info.ny; ++i) {
    for (int j = 0; j < info.nx; ++j) {

      const double dT         = T[info.idx2d(i, j)] - chpV.T_min;
      const double eps2_fit   = deps2_dT  * dT + chpV.eps2_min;
      const double sigma_fit  = dsigma_dT * dT + chpV.sigma_min;
      chpV.eps_2[info.idx2d(i, j)] = std::min( std::max( eps2_fit  , chpV.eps2_min )  , chpV.eps2_max );
      chpV.sigma[info.idx2d(i, j)] = std::min( std::max( sigma_fit , chpV.sigma_min ) , chpV.sigma_max );

    }
  }
  
  return chpV;

}
