#include <vector>
#include <algorithm>
#include "chparams.h"


double CHparamsScalar::compute_stability_limit(double dx , double dy) {
  double dmin = std::min( dx , dy );
  return 0.5 * dmin * dmin * dmin * dmin / eps_2;
};

double CHparamsVector::compute_stability_limit(double dx , double dy) {
  double dmin  = std::min( dx , dy );
  int idx_gmax = std::distance( eps_2.begin() , std::max_element( eps_2.begin() , eps_2.end() , abs_compare ) );
  double gmax  = eps_2[ idx_gmax ];
  return 0.5 * dmin * dmin * dmin * dmin / gmax;
};

double CHparamsVector::convert_temperature_to_flory_huggins( const double T ,
                                                             const double T_min ,
                                                             const double T_max ,
                                                             const double X_min ,
                                                             const double X_max ) {

  const double dX_dTinv   = ( X_max  - X_min ) / ( 1.0 / T_min - 1.0 / T_max );  
  const double dTinv      = 1.0 / T - 1.0 / T_max;
  const double X          = dX_dTinv * dTinv + X_min;

  return X;

};

double CHparamsVector::compute_eps2_from_polymer_params( const double X ,
                                                         const double m ,
                                                         const double L_kuhn ,
                                                         const double N ) {
  
  const double m_scaled   = 0.5 * ( 1.0 - m );
  const double Eps_2      = L_kuhn * L_kuhn / ( 3.0 * m_scaled * (1.0 - m_scaled) * (1.0 - m_scaled) * L_kuhn * L_kuhn * X * N * N );

  return Eps_2;
};

double CHparamsVector::compute_sigma_from_polymer_params( const double X ,
                                                          const double m ,
                                                          const double L_kuhn ,
                                                          const double L_omega ,
                                                          const double N  ) {
  
  const double m_scaled   = 0.5 * ( 1.0 - m );
  const double Sigma      = 36.0 * L_omega * L_omega / ( m_scaled * m_scaled * (1.0 - m_scaled) * (1.0 - m_scaled) * L_kuhn * L_kuhn * X * N * N );

  return Sigma;
};

void CHparamsVector::compute_and_set_eps2_and_sigma_from_polymer_params( const double T ,
                                                                         SimInfo& info ) {
  #pragma omp parallel for
  for (int i = 0; i < info.nx*info.ny; i++ ) {
    const double X = convert_temperature_to_flory_huggins( T , T_min , T_max , X_min , X_max );
    eps_2[i]       = compute_eps2_from_polymer_params( X , m[i] , L_kuhn , N );
    sigma[i]       = compute_sigma_from_polymer_params( X , m[i] , L_kuhn , L_omega , N );
  }

};

int SimInfo::idx2d_impl(int i, int j) {
  return j * ny + i;
};

// regular modulo operator gives negative values without this
int SimInfo::mod(int a, int b) {
  return (a % b + b) % b;
};

int SimInfo::idx2d(int i, int j)
{
  // modify the indices to map to a periodic mesh. need two levels for the 4th order operator.
  // i coordinates:
  i = mod(i, ny);
  j = mod(j, nx);
  
  return idx2d_impl(i, j);
};

