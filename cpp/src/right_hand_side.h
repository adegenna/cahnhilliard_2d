#ifndef __RIGHT_HAND_SIDE_H__
#define __RIGHT_HAND_SIDE_H__

#include <vector>

class RightHandSide {
public:
  virtual void rhs(const std::vector<double> &c, std::vector<double> &dcdt,
                   const double t) = 0;
  virtual void write_state( const std::vector<double> &x , const int idx , const int nx , const int ny ) = 0;
  virtual void setInitialConditions(std::vector<double> &x) = 0;
  void operator()(const std::vector<double> &c, std::vector<double> &dcdt, const double t)
  {
    rhs(c,dcdt,t);
  }
  double l2residual(const std::vector<double> &c);
  
};

#endif
