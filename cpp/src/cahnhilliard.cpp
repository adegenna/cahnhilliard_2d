#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <omp.h>
#include <boost/numeric/odeint.hpp>
#include "cahnhilliard.h"

CahnHilliard2DRHS::CahnHilliard2DRHS(CHparams& chp)
  : D_(chp.m), gamma_(chp.gam), b_(chp.b), u_(chp.u), alpha_(chp.alpha), phi_star_(chp.phi_star), nx_(chp.nx), dx_(chp.dx), sigma_(chp.sigma) , noise_dist_(0.0,1.0) , param_type_(chp.param_type), D_xy_(chp.m_xy), gamma_xy_(chp.gam_xy), b_xy_(chp.b_xy), u_xy_(chp.u_xy), alpha_xy_(chp.alpha_xy), phi_star_xy_(chp.phi_star_xy), sigma_xy_(chp.sigma_xy)
  {
    std::cout << "Initialized Cahn-Hilliard equation with ";
    if (param_type_ == 0)
      std::cout << "scalar parameters" << std::endl;
    else
      std::cout << "spatial-field parameters" << std::endl;
 
  }

  /*
  Cahn-Hilliard:
  
  dc/dt = D*laplacian( u*c^3 - b*c) - D*gamma*biharm(c)
  
  expanding out RHS into individual differentials:
  D*laplacian( u*c^3 - b*c) - D*gamma*biharm(c)
  assuming constant gamma.

  need a d^4 and a d^2 operator.
  */

CahnHilliard2DRHS::~CahnHilliard2DRHS() { };

void CahnHilliard2DRHS::rhs_scalar_parameters(const state_type &c, state_type &dcdt, const double t)
  {
    dcdt.resize(nx_*nx_);

    // evaluate the second order term, 5 point central stencil
    # pragma omp parallel for
    for (int i = 0; i < nx_; ++i)
    {
      for (int j = 0; j < nx_; ++j)
      {
        const double c_i   = laplace_component(c[idx2d(i, j)]);
        const double c_im1 = laplace_component(c[idx2d(i - 1, j)]);
        const double c_ip1 = laplace_component(c[idx2d(i + 1, j)]);
        const double c_jm1 = laplace_component(c[idx2d(i, j - 1)]);
        const double c_jp1 = laplace_component(c[idx2d(i, j + 1)]);
        dcdt[idx2d(i, j)]  = (D_ / (dx_ * dx_)) * (c_im1 + c_ip1 + c_jm1 + c_jp1 - 4.0 * c_i);
      }
    }

    // evaluate the 4th order term, 9 point central stencil
    # pragma omp parallel for
    for (int i = 0; i < nx_; ++i){
      for (int j = 0; j < nx_; ++j){
        const double c_i   = c[idx2d(i, j)];
        const double c_im1 = c[idx2d(i - 1, j)];
        const double c_ip1 = c[idx2d(i + 1, j)];
        const double c_im2 = c[idx2d(i - 2, j)];
        const double c_ip2 = c[idx2d(i + 2, j)];
        const double c_jm1 = c[idx2d(i, j - 1)];
        const double c_jp1 = c[idx2d(i, j + 1)];
        const double c_jm2 = c[idx2d(i, j - 2)];
        const double c_jp2 = c[idx2d(i, j + 2)];
	const double c_ul  = c[idx2d(i-1 , j-1)];
	const double c_ur  = c[idx2d(i-1 , j+1)];
	const double c_bl  = c[idx2d(i+1 , j-1)];
	const double c_br  = c[idx2d(i+1 , j+1)];

        // x-direction u_xxxx
        dcdt[idx2d(i,j)] -= D_ * gamma_ /(dx_*dx_*dx_*dx_) * 
          (c_ip2 - 4.0*c_ip1 + 6.0*c_i - 4.0*c_im1 + c_im2);

        // y-direction u_yyyy
        dcdt[idx2d(i,j)] -= D_ * gamma_ /(dx_*dx_*dx_*dx_) * 
          (c_jp2 - 4.0*c_jp1 + 6.0*c_i - 4.0*c_jm1 + c_jm2);

	// mixed term 2*u_xxyy
	dcdt[idx2d(i,j)] -= D_ * gamma_ /(dx_*dx_*dx_*dx_) * 
          2 * (4*c_i - 2*(c_im1 + c_ip1 + c_jm1 + c_jp1) + c_ul + c_ur + c_bl + c_br );
      }
    }

    // evaluate linear term
    # pragma omp parallel for
    for (int i = 0; i < nx_; ++i){
      for (int j = 0; j < nx_; ++j){
        const double c_i   = c[idx2d(i, j)];
    	dcdt[idx2d(i,j)]  -= alpha_ * ( c_i - phi_star_ );
      }
    }

    // evaluate the noise term
    state_type noise(nx_*nx_);
    # pragma omp parallel for
    for (int i = 0; i < nx_; ++i){
      for (int j = 0; j < nx_; ++j){
        noise[idx2d(i,j)]  = noise_dist_(generator_);
      }
    }
    double mean_noise = std::accumulate( noise.begin() , noise.end() , 0.0 ) / noise.size();

    # pragma omp parallel for
    for (int i = 0; i < nx_; ++i){
      for (int j = 0; j < nx_; ++j){
        const double c_i   = c[idx2d(i, j)];
	dcdt[idx2d(i,j)]  += sigma_ * (noise[i,j] - mean_noise);
      }
    }
    
  }

void CahnHilliard2DRHS::rhs_field_parameters(const state_type &c, state_type &dcdt, const double t)
  {
    dcdt.resize(nx_*nx_);

    // evaluate the second order term, 5 point central stencil
    # pragma omp parallel for
    for (int i = 0; i < nx_; ++i)
    {
      for (int j = 0; j < nx_; ++j)
      {
        const double c_i   = laplace_component_field(c , i     , j);
        const double c_im1 = laplace_component_field(c , i - 1 , j);
        const double c_ip1 = laplace_component_field(c , i + 1 , j);
        const double c_jm1 = laplace_component_field(c , i     , j - 1);
        const double c_jp1 = laplace_component_field(c , i     , j + 1);
        dcdt[idx2d(i, j)]  = (D_xy_[idx2d(i, j)] / (dx_ * dx_)) * (c_im1 + c_ip1 + c_jm1 + c_jp1 - 4.0 * c_i);
      }
    }

    // evaluate the 4th order term, 9 point central stencil
    # pragma omp parallel for
    for (int i = 0; i < nx_; ++i){
      for (int j = 0; j < nx_; ++j){
        const double c_i   = c[idx2d(i, j)];
        const double c_im1 = c[idx2d(i - 1, j)];
        const double c_ip1 = c[idx2d(i + 1, j)];
        const double c_im2 = c[idx2d(i - 2, j)];
        const double c_ip2 = c[idx2d(i + 2, j)];
        const double c_jm1 = c[idx2d(i, j - 1)];
        const double c_jp1 = c[idx2d(i, j + 1)];
        const double c_jm2 = c[idx2d(i, j - 2)];
        const double c_jp2 = c[idx2d(i, j + 2)];
	const double c_ul  = c[idx2d(i-1 , j-1)];
	const double c_ur  = c[idx2d(i-1 , j+1)];
	const double c_bl  = c[idx2d(i+1 , j-1)];
	const double c_br  = c[idx2d(i+1 , j+1)];

        // x-direction u_xxxx
        dcdt[idx2d(i,j)] -= D_xy_[idx2d(i, j)] * gamma_xy_[idx2d(i, j)] /(dx_*dx_*dx_*dx_) * 
          (c_ip2 - 4.0*c_ip1 + 6.0*c_i - 4.0*c_im1 + c_im2);

        // y-direction u_yyyy
        dcdt[idx2d(i,j)] -= D_xy_[idx2d(i, j)] * gamma_xy_[idx2d(i, j)] /(dx_*dx_*dx_*dx_) * 
          (c_jp2 - 4.0*c_jp1 + 6.0*c_i - 4.0*c_jm1 + c_jm2);

	// mixed term 2*u_xxyy
	dcdt[idx2d(i,j)] -= D_xy_[idx2d(i, j)] * gamma_xy_[idx2d(i, j)] /(dx_*dx_*dx_*dx_) * 
          2 * (4*c_i - 2*(c_im1 + c_ip1 + c_jm1 + c_jp1) + c_ul + c_ur + c_bl + c_br );
      }
    }

    // evaluate linear term
    # pragma omp parallel for
    for (int i = 0; i < nx_; ++i){
      for (int j = 0; j < nx_; ++j){
        const double c_i   = c[idx2d(i, j)];
    	dcdt[idx2d(i,j)]  -= alpha_xy_[idx2d(i, j)] * ( c_i - phi_star_xy_[idx2d(i, j)] );
      }
    }

    // evaluate the noise term
    state_type noise(nx_*nx_);
    # pragma omp parallel for
    for (int i = 0; i < nx_; ++i){
      for (int j = 0; j < nx_; ++j){
        noise[idx2d(i,j)]  = noise_dist_(generator_);
      }
    }
    double mean_noise = std::accumulate( noise.begin() , noise.end() , 0.0 ) / noise.size();

    # pragma omp parallel for
    for (int i = 0; i < nx_; ++i){
      for (int j = 0; j < nx_; ++j){
        const double c_i   = c[idx2d(i, j)];
	dcdt[idx2d(i,j)]  += sigma_xy_[idx2d(i, j)] * (noise[i,j] - mean_noise);
      }
    }
    
  }

void CahnHilliard2DRHS::operator()(const state_type &c, state_type &dcdt, const double t)
{
  if (param_type_ == 0)
    rhs_scalar_parameters(c,dcdt,t);
  else
    rhs_field_parameters(c,dcdt,t);
}

void CahnHilliard2DRHS::setInitialConditions(state_type &x)
  {
    x.resize(nx_ * nx_);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-1.0,1.0);

    // double initial_value = -1.0;
    for (int i = 0; i < nx_; ++i)
    {
      for (int j = 0; j < nx_; ++j)
      {
        x[idx2d(i,j)] = distribution(generator) * 0.005;
      }
    }
  }

double CahnHilliard2DRHS::l2residual(const state_type&c)
  {
    state_type dcdt;
    (*this)(c, dcdt, 0);
    double res = 0;
    for (int i = 0; i < nx_*nx_; ++i){
      res += dcdt[i] * dcdt[i];
    }
    return sqrt(res);
  }

double CahnHilliard2DRHS::laplace_component(double c)
  {
    return D_ * u_ * (c * c * c) - D_ * b_ * c;
  }

double CahnHilliard2DRHS::laplace_component_field(const state_type& c, int i, int j)
  {
    return D_xy_[idx2d(i, j)] * u_xy_[idx2d(i, j)] *
      (c[idx2d(i, j)] * c[idx2d(i, j)] * c[idx2d(i, j)]) -
      D_xy_[idx2d(i, j)] * b_xy_[idx2d(i, j)] * c[idx2d(i, j)];
  }

int CahnHilliard2DRHS::idx2d_impl(int i, int j)
  {
    return i * nx_ + j;
  }
  
  // regular modulo operator gives negative values without this
int CahnHilliard2DRHS::mod(int a, int b)
  { return (a%b+b)%b; }

int CahnHilliard2DRHS::idx2d(int i, int j)
  {
    // modify the indices to map to a periodic mesh. need two levels for the 4th order operator.
    // i coordinates:
    i = mod(i, nx_);
    j = mod(j, nx_);

    return idx2d_impl(i, j);
  }

struct Recorder
{
  int nx;
  std::ofstream& out;
  
  Recorder( std::ofstream& out , int nx )
    : out( out ) , nx( nx ) { }

  void operator()( const state_type &x , const double t )
  {
    for (int i = 0; i < nx; ++i){
      for (int j = 0; j < nx; ++j){
	out << x[i * nx + j] << " ";
      }
      out << std::endl;
    }
  }
  
};


void write_state(const state_type &x , const int idx , const int nx )
{
  std::ofstream out;
  out.open( "C_" + std::to_string(idx) + ".out" );
  out.precision(16);
  
  for (int i = 0; i < nx; ++i){
    for (int j = 0; j < nx; ++j){
      out << x[i * nx + j] << " ";
    }
  }

  out.close();
};

void run_ch_solver_checkpointing(CHparams& chparams)
{
  CahnHilliard2DRHS rhs(chparams);

  state_type x;
  rhs.setInitialConditions(x);

  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<state_type> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  double time                  = 0.0;
  const double stability_limit = 0.5*chparams.dx*chparams.dx*chparams.dx*chparams.dx/chparams.m/chparams.gam; // just an estimate
  double dt_initial            = stability_limit * 0.5;
  
  const double res0            = rhs.l2residual(x);

  std::cout << "residual at initial condition: " << res0 << std::endl;
  write_state(x,0,chparams.nx);

  double t_steps = int( chparams.tf / chparams.dt_check );
  for (int i = 0; i < t_steps; ++i){
    integrate_adaptive(controlled_stepper, rhs, x, time, time + chparams.dt_check, dt_initial);
    time += chparams.dt_check;
    std::cout << "iter: " << i+1 << " , t = " << time << ", relative residual: " << rhs.l2residual(x) / res0 << std::endl;
    write_state(x,i+1,chparams.nx);
  }

};

void run_ch_solver(CHparams& chparams)
{
  CahnHilliard2DRHS rhs(chparams);

  state_type x;
  if (chparams.t0 == 0) {
    rhs.setInitialConditions(x);
    int iter = 0;
  }
  else {
    x        = chparams.x;
    int iter = chparams.iter;
  }

  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<state_type> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  const double stability_limit = 0.5*chparams.dx*chparams.dx*chparams.dx*chparams.dx/chparams.m/chparams.gam; // just an estimate
  double dt_initial            = stability_limit * 0.5;
  double dt_check_residual     = dt_initial * 10.0;
  
  const double res0            = rhs.l2residual(x);

  std::cout << "residual at initial condition: " << res0 << std::endl;
  if (chparams.iter == 0)
    write_state(x,0,chparams.nx);

  integrate_adaptive(controlled_stepper, rhs, x, chparams.t0, chparams.tf, dt_initial);
  
  chparams.iter += 1;
  std::cout << "iter: " << chparams.iter << " , t = " << chparams.tf << ", relative residual: " << rhs.l2residual(x) / res0 << std::endl;
  write_state(x,chparams.iter,chparams.nx);
  chparams.x = x;

};
