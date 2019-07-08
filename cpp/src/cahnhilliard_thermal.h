#ifndef __CAHNHILLIARD_THERMAL_H__
#define __CAHNHILLIARD_THERMAL_H__

#include <vector>
#include <algorithm>
#include "chparams.h"

class CahnHilliard2DRHS_thermal {

 public:

  CahnHilliard2DRHS_thermal(CHparamsScalar& chp , SimInfo& info);
  CahnHilliard2DRHS_thermal(CHparamsVector& chp , SimInfo& info);
  ~CahnHilliard2DRHS_thermal();
  void rhs(const std::vector<double> &c, std::vector<double> &dcdt, const double t);
  void operator()(const std::vector<double> &c, std::vector<double> &dcdt, const double t);
  void setInitialConditions(std::vector<double> &x);
  double l2residual(const std::vector<double> &c);
  void write_state( const std::vector<double> &x , const int idx , const int nx , const int ny );
  
 private:

  CHparamsVector chpV_;
  SimInfo& info_;
  void (*ch_rhs_) (const std::vector<double>&, std::vector<double>&, double, CHparamsVector&, SimInfo&);

  std::default_random_engine generator_;
  std::normal_distribution<double> noise_dist_;
    
};





#endif
