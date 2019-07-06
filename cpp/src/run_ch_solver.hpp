#include <boost/numeric/odeint.hpp>
#include "stochastic_euler.hpp"
#include "cahnhilliard.h"
#include "cahnhilliard_thermal.h"

template<typename T_params , typename T_rhs>
void run_ch_solver( T_params& chparams , SimInfo& info , T_rhs& rhs )
{
  
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

  const double stability_limit = chparams.compute_stability_limit(info.dx , info.dy); // just an estimate
  const double res0            = rhs.l2residual(x);

  std::cout << "residual at initial condition: " << res0 << std::endl;
  if (info.iter == 0)
    rhs.write_state(x,0,info.nx,info.ny);

  if (chparams.sigma_noise < 1e-2) {
    std::cout << "Solving deterministic (noise-free) CH" << std::endl;
    integrate_adaptive(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
    //boost::numeric::odeint::integrate_const(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
  }
  else {
    std::cout << "Solving stochastic CH" << std::endl;
    boost::mt19937 rng;
    boost::numeric::odeint::integrate_const( stochastic_euler() ,
		     std::make_pair( rhs , ornstein_stoch( rng , chparams.sigma_noise ) ),
		     x , info.t0 , info.tf , stability_limit/40. );
  }
  info.iter += 1;
  std::cout << "iter: " << info.iter << " , t = " << info.tf << ", relative residual: " << rhs.l2residual(x) / res0 << std::endl;
  rhs.write_state(x,info.iter,info.nx,info.ny);
  info.x = x;

};
