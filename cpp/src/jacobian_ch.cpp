#include "jacobian_ch.h"

std::vector<double> evaluate_Jv_of_deterministic_ch( const std::vector<double>& v ,
						     CHparamsVector& chpV , 
						     SimInfo& info ) {

  // Computes action of Jacobian(F) on a vector for deterministic CH

  std::vector<double> Jv( v.size() );

  // Term -eps^2 * \nabla^4 (v)
  # pragma omp parallel for
  for (int i = 0; i < info.ny; ++i) {
    for (int j = 0; j < info.nx; ++j) {
      
      const double v_i   = v[info.idx2d(i, j)];
      const double v_im1 = v[info.idx2d(i - 1, j)];
      const double v_ip1 = v[info.idx2d(i + 1, j)];
      const double v_im2 = v[info.idx2d(i - 2, j)];
      const double v_ip2 = v[info.idx2d(i + 2, j)];
      const double v_jm1 = v[info.idx2d(i, j - 1)];
      const double v_jp1 = v[info.idx2d(i, j + 1)];
      const double v_jm2 = v[info.idx2d(i, j - 2)];
      const double v_jp2 = v[info.idx2d(i, j + 2)];
      const double v_ul  = v[info.idx2d(i-1 , j-1)];
      const double v_ur  = v[info.idx2d(i-1 , j+1)];
      const double v_bl  = v[info.idx2d(i+1 , j-1)];
      const double v_br  = v[info.idx2d(i+1 , j+1)];

      // y-direction u_yyyy
      double dyyyy = 1.0 / (info.dy * info.dy * info.dy * info.dy) * 
        (v_ip2 - 4.0*v_ip1 + 6.0*v_i - 4.0*v_im1 + v_im2);

      // x-direction u_xxxx
      double dxxxx = 1.0 / (info.dx * info.dx * info.dx * info.dx) * 
        (v_jp2 - 4.0*v_jp1 + 6.0*v_i - 4.0*v_jm1 + v_jm2);

      // mixed term 2*u_xxyy
      double dxxyy = 1.0 / (info.dx * info.dx * info.dy * info.dy) * 
        2 * (4*v_i - 2*(v_im1 + v_ip1 + v_jm1 + v_jp1) + v_ul + v_ur + v_bl + v_br );

      Jv[info.idx2d(i,j)] = -chpV.eps_2[info.idx2d(i,j)] * ( dxxxx + dyyyy + dxxyy );
      
    }
  }

  // Term \nabla^2( u * 3v^2 - b * v )
  # pragma omp parallel for
  for (int i = 0; i < info.ny; ++i) {
    for (int j = 0; j < info.nx; ++j) {
      
      const double v_i   = v[info.idx2d(i, j)];
      const double v_im1 = v[info.idx2d(i - 1, j)];
      const double v_ip1 = v[info.idx2d(i + 1, j)];
      const double v_jm1 = v[info.idx2d(i, j - 1)];
      const double v_jp1 = v[info.idx2d(i, j + 1)];
      
      const double u_i   = chpV.u[info.idx2d(i, j)];
      const double u_im1 = chpV.u[info.idx2d(i - 1, j)];
      const double u_ip1 = chpV.u[info.idx2d(i + 1, j)];
      const double u_jm1 = chpV.u[info.idx2d(i, j - 1)];
      const double u_jp1 = chpV.u[info.idx2d(i, j + 1)];
      
      const double b_i   = chpV.b[info.idx2d(i, j)];
      const double b_im1 = chpV.b[info.idx2d(i - 1, j)];
      const double b_ip1 = chpV.b[info.idx2d(i + 1, j)];
      const double b_jm1 = chpV.b[info.idx2d(i, j - 1)];
      const double b_jp1 = chpV.b[info.idx2d(i, j + 1)];
      
      double dxx = (1.0 / (info.dx * info.dx)) * (
	3 * ( u_jp1*v_jp1*v_jp1 + u_jm1*v_jm1*v_jm1 - 2.0*u_i*v_i*v_i ) + 
	-   ( b_jp1*v_jp1       + b_jm1*v_jm1       - 2.0*b_i*v_i     ) );
      double dyy = (1.0 / (info.dy * info.dy)) * (
	3 * ( u_ip1*v_ip1*v_ip1 + u_im1*v_im1*v_im1 - 2.0*u_i*v_i*v_i )
	-   ( b_ip1*v_ip1       + b_im1*v_im1       - 2.0*b_i*v_i     ) );
      Jv[info.idx2d(i, j)] += dxx + dyy;
      
    }
  }

  // Term -\sigma * v
  # pragma omp parallel for
  for (int i = 0; i < info.ny; ++i){
    for (int j = 0; j < info.nx; ++j){
      
      const double v_i      = v[info.idx2d(i, j)];
      Jv[info.idx2d(i,j)]  += -chpV.sigma[info.idx2d(i,j)] * c_i;
      
    }
  }
  
  return Jv;
}
