#ifndef __CAHNHILLIARD_H__
#define __CAHNHILLIARD_H__

#include <vector>
#include <algorithm>

// replaced with lambda
// static bool abs_compare(int a, int b)
// {
//   return (std::abs(a) < std::abs(b));
// }

class SimInfo
{
  
 public:

  // these do nothing, and the compiler auto-generates these
  // SimInfo() { };
  // ~SimInfo() { };
  
  double t0;
  double tf;
  int iter = 0;
  double dx;
  int nx;
  double dt_check;
  std::vector<double> x;

};

class CHparamsScalar
{
  
 public:

  // these do nothing, and the compiler auto-generates these
  // CHparamsScalar()  { };
  // ~CHparamsScalar() { };

  double D;
  double gamma;
  double b;
  double u;
  double alpha;
  double phi_star;
  double sigma;

  double compute_stability_limit(double dx)
  {
    return 0.5 * dx * dx * dx * dx / D / gamma;
  };
  
};

class CHparamsVector
{
  
 public:

  // these do nothing, and the compiler auto-generates these
  // CHparamsVector()  { };
  // ~CHparamsVector() { };
  
  std::vector<double> D;
  std::vector<double> gamma;
  std::vector<double> b;
  std::vector<double> u;
  std::vector<double> alpha;
  std::vector<double> phi_star;
  double sigma;

  double compute_stability_limit(double dx)
  {
    // int idx_dmax = std::distance( D.begin()     , std::max_element( D.begin()     , D.end()     , abs_compare ) );
    // int idx_gmax = std::distance( gamma.begin() , std::max_element( gamma.begin() , gamma.end() , abs_compare ) );
    // double dmax = D[     idx_dmax ];
    // double gmax = gamma[ idx_gmax ];

    // lambda one-liner
    const double dmax = *std::max_element(D.begin(), D.end(), [](const double& a, const double& b){ return std::abs(a) < std::abs(b); } );
    const double gmax = *std::max_element(gamma.begin(), gamma.end(), [](const double& a, const double& b){ return std::abs(a) < std::abs(b); } );
    return 0.5 * dx * dx * dx * dx / dmax / gmax;
  };

  
};

class CahnHilliard2DRHS {

 public:

  CahnHilliard2DRHS(SimInfo& info);
  ~CahnHilliard2DRHS();


  void operator()(const std::vector<double> &c, std::vector<double> &dcdt, const double t);

  void setInitialConditions(std::vector<double> &x); // never used, remove? gets set by the integrator anyway

  double l2residual(const std::vector<double> &c);
  double get_gamma() const { return get_gamma_impl(); };
  
 protected:

  virtual void get_ij_values(int i , int j , CHparamsScalar& ch_ij) = 0;

  int idx2d(int i, int j);

 private:

// the point of having the derived classes is so that CahnHilliard2DRHS doesn't know or care about whether
// it's using vector or scalar parameters. it shouldn't own these, they should stay in their own classes.
// you want to depend on the abstraction, not the actual thing.
  // CHparamsScalar chpS_;
  // CHparamsVector chpV_;
  SimInfo info_;

  std::default_random_engine generator_;
  std::normal_distribution<double> noise_dist_;

  void rhs(const std::vector<double> &c, std::vector<double> &dcdt, const double t);
  virtual double get_gamma_impl() const = 0;
  
  double laplace_component(double c , double D , double u , double b);
  int idx2d_impl(int i, int j);
  int mod(int a, int b);
  CHparamsScalar ch_ij_;
  
};

class CahnHilliard2DRHS_Scalar : public CahnHilliard2DRHS {

 public:
  
  CahnHilliard2DRHS_Scalar(CHparamsScalar& chp , SimInfo& info)
    : CahnHilliard2DRHS(info),
      chpS_(chp)
      {}

 private:

  CHparamsScalar chpS_;

  void get_ij_values(int i , int j , CHparamsScalar& ch_ij) override;
  double get_gamma_impl() const override { return chpS_.gamma; }

};

class CahnHilliard2DRHS_Vector : public CahnHilliard2DRHS {

 public:
  
  CahnHilliard2DRHS_Vector(CHparamsVector& chp , SimInfo& info)
    : CahnHilliard2DRHS(info),
       chpV_(chp)
      {}

 private:

  CHparamsVector chpV_;

  void get_ij_values(int i , int j , CHparamsScalar& ch_ij) override;
  double get_gamma_impl() const override { return chpV_.gamma; }

};


void write_state(const std::vector<double> &x , const int idx , const int nx );

//void run_ch_solver_checkpointing(CHparams& chparams);




#endif
