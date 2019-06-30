#include <vector>
#include <algorithm>
#include "chparams.h"


double CHparamsScalar::compute_stability_limit(double dx) {
  return 0.5 * dx * dx * dx * dx / eps_2;
};

double CHparamsVector::compute_stability_limit(double dx) {
  int idx_gmax = std::distance( eps_2.begin() , std::max_element( eps_2.begin() , eps_2.end() , abs_compare ) );
  double gmax = eps_2[ idx_gmax ];
  return 0.5 * dx * dx * dx * dx / gmax;
};
