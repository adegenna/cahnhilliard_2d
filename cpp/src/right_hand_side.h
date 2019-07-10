#ifndef __RIGHT_HAND_SIDE_H__
#define __RIGHT_HAND_SIDE_H__

#include <vector>

class RightHandSide {
public:
  virtual void rhs(const std::vector<double> &c, std::vector<double> &dcdt,
                   const double t) = 0;


  void operator()(const std::vector<double> &c, std::vector<double> &dcdt, const double t)
  {
    rhs(c,dcdt,t);
  }

  double l2residual(const std::vector<double> &c);
};

#endif
