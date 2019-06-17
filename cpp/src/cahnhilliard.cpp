#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <omp.h>
#include <boost/numeric/odeint.hpp>
#include "cahnhilliard.h"
#include "stochastic_euler.hpp"


  /*
  Cahn-Hilliard:
  
  dc/dt = D*laplacian( u*c^3 - b*c) - D*gamma*biharm(c)
  
  expanding out RHS into individual differentials:
  D*laplacian( u*c^3 - b*c) - D*gamma*biharm(c)
  assuming constant gamma.

  need a d^4 and a d^2 operator.
  */


CahnHilliard2DRHS::CahnHilliard2DRHS(CHparams& chp)
  : noise_dist_(0.0,1.0) , chp_(chp)
  {
    std::cout << "Initialized Cahn-Hilliard equation with ";
    if (chp_.param_type == 0)
      std::cout << "scalar parameters" << std::endl;
    else
      std::cout << "spatial-field parameters" << std::endl;
 
  }

CahnHilliard2DRHS::~CahnHilliard2DRHS() { };

void CahnHilliard2DRHS::rhs_scalar_parameters(const std::vector<double> &c, std::vector<double> &dcdt, const double t)
  {
    dcdt.resize(chp_.nx*chp_.nx);

    // evaluate the second order term, 5 point central stencil
    # pragma omp parallel for
    for (int i = 0; i < chp_.nx; ++i)
    {
      for (int j = 0; j < chp_.nx; ++j)
      {
        const double c_i   = laplace_component(c[idx2d(i, j)]);
        const double c_im1 = laplace_component(c[idx2d(i - 1, j)]);
        const double c_ip1 = laplace_component(c[idx2d(i + 1, j)]);
        const double c_jm1 = laplace_component(c[idx2d(i, j - 1)]);
        const double c_jp1 = laplace_component(c[idx2d(i, j + 1)]);
        dcdt[idx2d(i, j)]  = (chp_.D / (chp_.dx * chp_.dx)) * (c_im1 + c_ip1 + c_jm1 + c_jp1 - 4.0 * c_i);
      }
    }

    // evaluate the 4th order term, 9 point central stencil
    # pragma omp parallel for
    for (int i = 0; i < chp_.nx; ++i){
      for (int j = 0; j < chp_.nx; ++j){
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
        dcdt[idx2d(i,j)] -= chp_.D * chp_.gamma /(chp_.dx*chp_.dx*chp_.dx*chp_.dx) * 
          (c_ip2 - 4.0*c_ip1 + 6.0*c_i - 4.0*c_im1 + c_im2);

        // y-direction u_yyyy
        dcdt[idx2d(i,j)] -= chp_.D * chp_.gamma /(chp_.dx*chp_.dx*chp_.dx*chp_.dx) * 
          (c_jp2 - 4.0*c_jp1 + 6.0*c_i - 4.0*c_jm1 + c_jm2);

        // mixed term 2*u_xxyy
        dcdt[idx2d(i,j)] -= chp_.D * chp_.gamma /(chp_.dx*chp_.dx*chp_.dx*chp_.dx) * 
          2 * (4*c_i - 2*(c_im1 + c_ip1 + c_jm1 + c_jp1) + c_ul + c_ur + c_bl + c_br );
      }
    }

    // evaluate linear term
    # pragma omp parallel for
    for (int i = 0; i < chp_.nx; ++i){
      for (int j = 0; j < chp_.nx; ++j){
        const double c_i   = c[idx2d(i, j)];
        dcdt[idx2d(i,j)]  -= chp_.alpha * ( c_i - chp_.phi_star );
      }
    }
    
  }

void CahnHilliard2DRHS::rhs_field_parameters(const std::vector<double> &c, std::vector<double> &dcdt, const double t)
  {
    dcdt.resize(chp_.nx*chp_.nx);

    // evaluate the second order term, 5 point central stencil
    # pragma omp parallel for
    for (int i = 0; i < chp_.nx; ++i)
    {
      for (int j = 0; j < chp_.nx; ++j)
      {
        const double c_i   = laplace_component_field(c , i     , j);
        const double c_im1 = laplace_component_field(c , i - 1 , j);
        const double c_ip1 = laplace_component_field(c , i + 1 , j);
        const double c_jm1 = laplace_component_field(c , i     , j - 1);
        const double c_jp1 = laplace_component_field(c , i     , j + 1);
        dcdt[idx2d(i, j)]  = (chp_.D_xy[idx2d(i, j)] / (chp_.dx * chp_.dx)) * (c_im1 + c_ip1 + c_jm1 + c_jp1 - 4.0 * c_i);
      }
    }

    // evaluate the 4th order term, 9 point central stencil
    # pragma omp parallel for
    for (int i = 0; i < chp_.nx; ++i){
      for (int j = 0; j < chp_.nx; ++j){
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
        dcdt[idx2d(i,j)] -= chp_.D_xy[idx2d(i, j)] * chp_.gamma_xy[idx2d(i, j)] /(chp_.dx*chp_.dx*chp_.dx*chp_.dx) * 
          (c_ip2 - 4.0*c_ip1 + 6.0*c_i - 4.0*c_im1 + c_im2);

        // y-direction u_yyyy
        dcdt[idx2d(i,j)] -= chp_.D_xy[idx2d(i, j)] * chp_.gamma_xy[idx2d(i, j)] /(chp_.dx*chp_.dx*chp_.dx*chp_.dx) * 
          (c_jp2 - 4.0*c_jp1 + 6.0*c_i - 4.0*c_jm1 + c_jm2);

        // mixed term 2*u_xxyy
        dcdt[idx2d(i,j)] -= chp_.D_xy[idx2d(i, j)] * chp_.gamma_xy[idx2d(i, j)] /(chp_.dx*chp_.dx*chp_.dx*chp_.dx) * 
          2 * (4*c_i - 2*(c_im1 + c_ip1 + c_jm1 + c_jp1) + c_ul + c_ur + c_bl + c_br );
      }
    }

    // evaluate linear term
    # pragma omp parallel for
    for (int i = 0; i < chp_.nx; ++i){
      for (int j = 0; j < chp_.nx; ++j){
        const double c_i   = c[idx2d(i, j)];
        dcdt[idx2d(i,j)]  -= chp_.alpha_xy[idx2d(i, j)] * ( c_i - chp_.phi_star_xy[idx2d(i, j)] );
      }
    }
    
  }

void CahnHilliard2DRHS::operator()(const std::vector<double> &c, std::vector<double> &dcdt, const double t)
{
  if (chp_.param_type == 0)
    rhs_scalar_parameters(c,dcdt,t);
  else
    rhs_field_parameters(c,dcdt,t);
}

void CahnHilliard2DRHS::setInitialConditions(std::vector<double> &x)
  {
    x.resize(chp_.nx * chp_.nx);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-1.0,1.0);

    // double initial_value = -1.0;
    for (int i = 0; i < chp_.nx; ++i)
    {
      for (int j = 0; j < chp_.nx; ++j)
      {
        x[idx2d(i,j)] = distribution(generator) * 0.005;
      }
    }
  }

double CahnHilliard2DRHS::l2residual(const std::vector<double>&c)
  {
    std::vector<double> dcdt;
    (*this)(c, dcdt, 0);
    double res = 0;
    for (int i = 0; i < chp_.nx*chp_.nx; ++i){
      res += dcdt[i] * dcdt[i];
    }
    return sqrt(res);
  }

double CahnHilliard2DRHS::laplace_component(double c)
  {
    return chp_.D * chp_.u * (c * c * c) - chp_.D * chp_.b * c;
  }

double CahnHilliard2DRHS::laplace_component_field(const std::vector<double>& c, int i, int j)
  {
    return chp_.D_xy[idx2d(i, j)] * chp_.u_xy[idx2d(i, j)] *
      (c[idx2d(i, j)] * c[idx2d(i, j)] * c[idx2d(i, j)]) -
      chp_.D_xy[idx2d(i, j)] * chp_.b_xy[idx2d(i, j)] * c[idx2d(i, j)];
  }

int CahnHilliard2DRHS::idx2d_impl(int i, int j)
  {
    return i * chp_.nx + j;
  }
  
  // regular modulo operator gives negative values without this
int CahnHilliard2DRHS::mod(int a, int b)
  { return (a%b+b)%b; }

int CahnHilliard2DRHS::idx2d(int i, int j)
  {
    // modify the indices to map to a periodic mesh. need two levels for the 4th order operator.
    // i coordinates:
    i = mod(i, chp_.nx);
    j = mod(j, chp_.nx);

    return idx2d_impl(i, j);
  }

struct Recorder
{
  int nx;
  std::ofstream& out;
  
  Recorder( std::ofstream& out , int nx )
    : out( out ) , nx( nx ) { }

  void operator()( const std::vector<double> &x , const double t )
  {
    for (int i = 0; i < nx; ++i){
      for (int j = 0; j < nx; ++j){
        out << x[i * nx + j] << " ";
      }
      out << std::endl;
    }
  }
  
};


void write_state(const std::vector<double> &x , const int idx , const int nx )
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

  std::vector<double> x;
  rhs.setInitialConditions(x);

  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  double time                  = 0.0;
  const double stability_limit = 0.5*chparams.dx*chparams.dx*chparams.dx*chparams.dx/chparams.D/chparams.gamma; // just an estimate
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

  std::vector<double> x;
  if (chparams.t0 == 0) {
    rhs.setInitialConditions(x);
    int iter = 0;
  }
  else {
    x        = chparams.x;
    int iter = chparams.iter;
  }

  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  const double stability_limit = 0.5*chparams.dx*chparams.dx*chparams.dx*chparams.dx/chparams.D/chparams.gamma; // just an estimate
  double dt_initial            = stability_limit * 0.5;
  double dt_check_residual     = dt_initial * 10.0;
  
  const double res0            = rhs.l2residual(x);

  std::cout << "residual at initial condition: " << res0 << std::endl;
  if (chparams.iter == 0)
    write_state(x,0,chparams.nx);

  if (chparams.sigma < 1e-2) {
    std::cout << "Solving deterministic (noise-free) CH" << std::endl;
    integrate_adaptive(controlled_stepper, rhs, x, chparams.t0, chparams.tf, stability_limit/2.);
  }
  else {
    std::cout << "Solving stochastic CH" << std::endl;
    boost::mt19937 rng;
    boost::numeric::odeint::integrate_const( stochastic_euler() ,
                                             std::make_pair( rhs , ornstein_stoch( rng , chparams.sigma ) ),
                                             x , chparams.t0 , chparams.tf , stability_limit/40. );
  }
  chparams.iter += 1;
  std::cout << "iter: " << chparams.iter << " , t = " << chparams.tf << ", relative residual: " << rhs.l2residual(x) / res0 << std::endl;
  write_state(x,chparams.iter,chparams.nx);
  chparams.x = x;

};
