#include "jacobian_ch.h"
#include "newton_krylov.h"
#include "utils_ch.h"

double l2_norm( std::vector<double>& x ) {

  double norm = 0.0;

  #pragma omp parallel for
  for (int i = 0; i < x.size(); i++)
    norm += x[i]*x[i];
  norm = pow(norm , 0.5);

  return norm;
  
};

std::vector<double> newton_krylov_deterministic_ch( const std::vector<double>& c  , 
						    std::vector<double>& x ,
						    const double epsilon ,
						    const int max_iter ,
						    CHparamsVector& chpV ,
						    SimInfo& info ) {

  std::vector<double> R  = compute_ch_residual( c , x , chpV , info );
  int counter = 1;

  while ( l2_norm( R ) > epsilon ) {
    v = gmres( evaluate_Jv_of_deterministic_ch , R , x , chpV , info );
    x = x - v;
    R = compute_ch_residual( c , x , chpV , info );
    if ( counter > max_iter ) {
      std::cout << "Error: newton-krylov did not converge\n";
      break;
    }
    counter += 1;
    
  }
  
  return R;

};
