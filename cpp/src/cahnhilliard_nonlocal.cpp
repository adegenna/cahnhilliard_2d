#include "cahnhilliard_nonlocal.h"

void compute_ch_nonlocal(const std::vector<double> &c,
			 std::vector<double> &dcdt,
			 const double t,
			 CHparamsVector chpV,
			 SimInfo& info) {

  // Computes deterministic nonlocal CH dynamics
  // dc/dt = laplacian( u*c^3 - b*c ) - eps_2*biharm(c) - sigma*(c - m)

  // evaluate the second order term, 5 point central stencil
  # pragma omp parallel for
  for (int i = 0; i < info.nx; ++i) {
    for (int j = 0; j < info.nx; ++j) {
        
      const double c_i   = laplace_component( info.idx2d(i, j)      , c , chpV.u , chpV.b );
      const double c_im1 = laplace_component( info.idx2d(i - 1, j)  , c , chpV.u , chpV.b );
      const double c_ip1 = laplace_component( info.idx2d(i + 1, j)  , c , chpV.u , chpV.b );
      const double c_jm1 = laplace_component( info.idx2d(i, j - 1)  , c , chpV.u , chpV.b );
      const double c_jp1 = laplace_component( info.idx2d(i, j + 1)  , c , chpV.u , chpV.b );
        
      dcdt[info.idx2d(i, j)]  = (1.0 / (info.dx * info.dx)) * (c_im1 + c_ip1 + c_jm1 + c_jp1 - 4.0 * c_i);
    }
  }

  // evaluate the 4th order term, 9 point central stencil
  # pragma omp parallel for
  for (int i = 0; i < info.nx; ++i) {
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

      // x-direction u_xxxx
      dcdt[info.idx2d(i,j)] -= chpV.eps_2[info.idx2d(i,j)] /(info.dx*info.dx*info.dx*info.dx) * 
        (c_ip2 - 4.0*c_ip1 + 6.0*c_i - 4.0*c_im1 + c_im2);

      // y-direction u_yyyy
      dcdt[info.idx2d(i,j)] -= chpV.eps_2[info.idx2d(i,j)] /(info.dx*info.dx*info.dx*info.dx) * 
        (c_jp2 - 4.0*c_jp1 + 6.0*c_i - 4.0*c_jm1 + c_jm2);

      // mixed term 2*u_xxyy
      dcdt[info.idx2d(i,j)] -= chpV.eps_2[info.idx2d(i,j)] /(info.dx*info.dx*info.dx*info.dx) * 
        2 * (4*c_i - 2*(c_im1 + c_ip1 + c_jm1 + c_jp1) + c_ul + c_ur + c_bl + c_br );
    }
  }

  // evaluate linear term
  # pragma omp parallel for
  for (int i = 0; i < info.nx; ++i){
    for (int j = 0; j < info.nx; ++j){
        
      const double c_i   = c[info.idx2d(i, j)];
        
      dcdt[info.idx2d(i,j)]  -= chpV.sigma[info.idx2d(i,j)] * ( c_i - chpV.m[info.idx2d(i,j)] );
    }
  }


}

double laplace_component(int i ,
                         const std::vector<double>& c ,
                         const std::vector<double>& u ,
                         const std::vector<double>& b ) {

  return u[i] * (c[i] * c[i] * c[i]) - b[i] * c[i];
}
