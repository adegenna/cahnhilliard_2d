#include "right_hand_side.h"

#include <cmath>

double RightHandSide::l2residual(const std::vector<double> &c) {
  std::vector<double> dcdt;
  (*this)(c, dcdt, 0);
  double res = 0;
  for (int i = 0; i < dcdt.size(); ++i) {
    res += dcdt[i] * dcdt[i];
  }
  return std::sqrt(res);
}
