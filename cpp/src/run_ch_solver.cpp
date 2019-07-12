#include <boost/numeric/odeint.hpp>
#include "stochastic_euler.hpp"
#include "cahnhilliard.h"
#include "cahnhilliard_thermal.h"
#include "cahnhilliard_thermal_nodiffusion.h"
#include "run_ch_solver.h"

void run_ch_solver_non_thermal( CHparamsVector& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS rhs = CahnHilliard2DRHS( chparams , info );
  
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


void run_ch_solver_thermal_no_diffusion( CHparamsVector& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS_thermal_nodiffusion rhs = CahnHilliard2DRHS_thermal_nodiffusion( chparams , info );
  
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


void run_ch_solver_thermal_with_diffusion( CHparamsVector& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS_thermal rhs = CahnHilliard2DRHS_thermal( chparams , info );
  
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


void run_ch_solver_non_thermal( CHparamsScalar& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS rhs = CahnHilliard2DRHS( chparams , info );
  
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


void run_ch_solver_thermal_no_diffusion( CHparamsScalar& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS_thermal_nodiffusion rhs = CahnHilliard2DRHS_thermal_nodiffusion( chparams , info );
  
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


void run_ch_solver_thermal_with_diffusion( CHparamsScalar& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS_thermal rhs = CahnHilliard2DRHS_thermal( chparams , info );
  
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

void run_ch_solver_vector( CHparamsVector& chparams , SimInfo& info )
{
  
  if      (info.rhs_type.compare("ch_non_thermal") == 0) {
    run_ch_solver_non_thermal( chparams , info );
  }
  else if (info.rhs_type.compare("ch_thermal_no_diffusion") == 0) {
    run_ch_solver_thermal_no_diffusion( chparams , info );
  }
  else if (info.rhs_type.compare("ch_thermal_with_diffusion") == 0) {
    run_ch_solver_thermal_with_diffusion( chparams , info );
  }

};

void run_ch_solver_scalar( CHparamsScalar& chparams , SimInfo& info )
{
  
  if      (info.rhs_type.compare("ch_non_thermal") == 0) {
    run_ch_solver_non_thermal( chparams , info );
  }
  else if (info.rhs_type.compare("ch_thermal_no_diffusion") == 0) {
    run_ch_solver_thermal_no_diffusion( chparams , info );
  }
  else if (info.rhs_type.compare("ch_thermal_with_diffusion") == 0) {
    run_ch_solver_thermal_with_diffusion( chparams , info );
  }

};
