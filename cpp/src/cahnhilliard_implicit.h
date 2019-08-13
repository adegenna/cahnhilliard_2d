#ifndef __CAHNHILLIARD_IMPLICIT_H__
#define __CAHNHILLIARD_IMPLICIT_H__

#include <vector>
#include <algorithm>
#include <random>
#include "chparams.h"
#include "right_hand_side.h"

class CahnHilliardImplicit2DRHS : public RightHandSide {

 public:

  CahnHilliardImplicit2DRHS(CHparamsScalar& chp , SimInfo& info , DM& da , TS& ts , Vec& F);
  CahnHilliardImplicit2DRHS(CHparamsVector& chp , SimInfo& info , DM& da , TS& ts , Vec& F);
  ~CahnHilliardImplicit2DRHS();
  void rhs(const std::vector<double> &c, std::vector<double> &dcdt, const double t) override;
  void write_state(const std::vector<double> &x , const int idx , const int nx , const int ny , std::string& outdir) override;
  void setInitialConditions(std::vector<double> &x);
  
 private:

  CHparamsVector chpV_;
  SimInfo& info_;
  void (*ch_rhs_) (const std::vector<double>&, std::vector<double>&, double, CHparamsVector&, SimInfo&, TS&, Vec&, DM&);

  std::default_random_engine generator_;
  std::normal_distribution<double> noise_dist_;

  // Petsc data members
  TS ts_;
  Vec F_;
  DM da_;
  
};

#endif
