#ifndef __CAHNHILLIARD_H__
#define __CAHNHILLIARD_H__

#include <vector>
#include <algorithm>
#include "chparams.h"
#include "right_hand_side.h"

class CahnHilliard2DRHS : public RightHandSide {

 public:

  CahnHilliard2DRHS(CHparamsScalar& chp , SimInfo& info);
  CahnHilliard2DRHS(CHparamsVector& chp , SimInfo& info);
  ~CahnHilliard2DRHS();
  void rhs(const std::vector<double> &c, std::vector<double> &dcdt, const double t) override;
  void setInitialConditions(std::vector<double> &x);
  void write_state(const std::vector<double> &x , const int idx , const int nx , const int ny );
  
 private:

  CHparamsVector chpV_;
  SimInfo& info_;
  void (*ch_rhs_) (const std::vector<double>&, std::vector<double>&, double, CHparamsVector&, SimInfo&);

  std::default_random_engine generator_;
  std::normal_distribution<double> noise_dist_;
    
};

#endif
