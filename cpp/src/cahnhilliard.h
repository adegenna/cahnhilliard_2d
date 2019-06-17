#ifndef __CAHNHILLIARD_H__
#define __CAHNHILLIARD_H__

#include <vector>
#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include "stochastic_euler.hpp"

static bool abs_compare(int a, int b)
{
  return (std::abs(a) < std::abs(b));
}

class SimInfo
{
  
 public:

  SimInfo() { };
  ~SimInfo() { };
  
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

  CHparamsScalar()  { };
  ~CHparamsScalar() { };

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

  CHparamsVector()  { };
  ~CHparamsVector() { };
  
  std::vector<double> D;
  std::vector<double> gamma;
  std::vector<double> b;
  std::vector<double> u;
  std::vector<double> alpha;
  std::vector<double> phi_star;
  double sigma;

  double compute_stability_limit(double dx)
  {
    int idx_dmax = std::distance( D.begin()     , std::max_element( D.begin()     , D.end()     , abs_compare ) );
    int idx_gmax = std::distance( gamma.begin() , std::max_element( gamma.begin() , gamma.end() , abs_compare ) );
    double dmax = D[     idx_dmax ];
    double gmax = gamma[ idx_gmax ];
    return 0.5 * dx * dx * dx * dx / dmax / gmax;
  };

  
};

class CahnHilliard2DRHS {

 public:

  CahnHilliard2DRHS(CHparamsScalar& chp , SimInfo& info);
  CahnHilliard2DRHS(CHparamsVector& chp , SimInfo& info);
  ~CahnHilliard2DRHS();
  void rhs(const std::vector<double> &c, std::vector<double> &dcdt, const double t);
  void operator()(const std::vector<double> &c, std::vector<double> &dcdt, const double t);
  void setInitialConditions(std::vector<double> &x);
  double l2residual(const std::vector<double> &c);
  
 protected:

  CHparamsScalar chpS_;
  CHparamsVector chpV_;
  SimInfo& info_;

  std::default_random_engine generator_;
  std::normal_distribution<double> noise_dist_;

  virtual void get_ij_values(int i , int j , CHparamsScalar& ch_ij) ;
  
  double laplace_component(double c , double D , double u , double b);
  int idx2d_impl(int i, int j);
  int mod(int a, int b);
  int idx2d(int i, int j);
  CHparamsScalar ch_ij_;
  
};

class CahnHilliard2DRHS_Scalar : public CahnHilliard2DRHS {

 public:
  CahnHilliard2DRHS_Scalar(CHparamsScalar& chp , SimInfo& info)
    : CahnHilliard2DRHS(chp , info)
      {}

 private:

  void get_ij_values(int i , int j , CHparamsScalar& ch_ij) ;

};

class CahnHilliard2DRHS_Vector : public CahnHilliard2DRHS {

 private:

  void get_ij_values(int i , int j , CHparamsScalar& ch_ij) ;

};


void write_state(const std::vector<double> &x , const int idx , const int nx );

//void run_ch_solver_checkpointing(CHparams& chparams);

template<typename T>
void run_ch_solver(T& chparams , SimInfo& info )
{

  auto rhs = CahnHilliard2DRHS(chparams,info);

  std::vector<double> x;
  if (info.t0 == 0) {
    rhs.setInitialConditions(x);
    int iter = 0;
  }
  else {
    x        = info.x;
    int iter = info.iter;
  }

  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  const double stability_limit = chparams.compute_stability_limit(info.dx); // just an estimate
  double dt_initial            = stability_limit * 0.5;
  double dt_check_residual     = dt_initial * 10.0;
  
  const double res0            = rhs.l2residual(x);

  std::cout << "residual at initial condition: " << res0 << std::endl;
  if (info.iter == 0)
    write_state(x,0,info.nx);

  if (chparams.sigma < 1e-2) {
    std::cout << "Solving deterministic (noise-free) CH" << std::endl;
    integrate_adaptive(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
  }
  else {
    std::cout << "Solving stochastic CH" << std::endl;
    boost::mt19937 rng;
    boost::numeric::odeint::integrate_const( stochastic_euler() ,
                                             std::make_pair( rhs , ornstein_stoch( rng , chparams.sigma ) ),
                                             x , info.t0 , info.tf , stability_limit/40. );
  }
  info.iter += 1;
  std::cout << "iter: " << info.iter << " , t = " << info.tf << ", relative residual: " << rhs.l2residual(x) / res0 << std::endl;
  write_state(x,info.iter,info.nx);
  info.x = x;

};



#endif
