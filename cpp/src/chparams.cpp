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

