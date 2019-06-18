#include <boost/numeric/odeint.hpp>
#include "stochastic_euler.hpp"
#include "cahnhilliard.h"

template<typename T>
void run_ch_solver(T& chparams , SimInfo& info )
{

  CahnHilliard2DRHS rhs = CahnHilliard2DRHS(chparams , info);

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
