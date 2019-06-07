#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <omp.h>

#include <boost/numeric/odeint.hpp>

//#include "matplotlibcpp.h"

typedef std::vector<double> state_type;

struct CHparams
{
  double m;
  double gam;
  double b;
  double u;
  double sig;
  double phi_star;
};


class CahnHilliard2DRHS
{
public:
  CahnHilliard2DRHS(CHparams& chp, int nx, double dx)
    : D_(chp.m), gamma_(chp.gam), b_(chp.b), u_(chp.u), sig_(chp.sig), phi_star_(chp.phi_star), nx_(nx), dx_(dx)
  {
    std::cout << "Initialized Cahn-Hilliard equation with D_ " << D_ 
      << " gamma_ " << gamma_ << " dx_ " << nx_ << " dx_ " << dx_ << std::endl;
  }

  /*
  Cahn-Hilliard:
  
  dc/dt = D*laplacian( u*c^3 - b*c) - D*gamma*biharm(c)
  
  expanding out RHS into individual differentials:
  D*laplacian( u*c^3 - b*c) - D*gamma*biharm(c)
  assuming constant gamma.

  need a d^4 and a d^2 operator.
  */
  void operator()(const state_type &c, state_type &dcdt, const double t)
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

        // x-direction
        dcdt[idx2d(i,j)] -= D_ * gamma_ /(dx_*dx_*dx_*dx_) * 
          (c_ip2 - 4.0*c_ip1 + 6.0*c_i - 4.0*c_im1 + c_im2);

        // y-direction
        dcdt[idx2d(i,j)] -= D_ * gamma_ /(dx_*dx_*dx_*dx_) * 
          (c_jp2 - 4.0*c_jp1 + 6.0*c_i - 4.0*c_jm1 + c_jm2);
      }
    }

    // evaluate linear term
    # pragma omp parallel for
    for (int i = 0; i < nx_; ++i){
      for (int j = 0; j < nx_; ++j){
        const double c_i   = c[idx2d(i, j)];
    	dcdt[idx2d(i,j)]  += sig_ * ( c_i - phi_star_ );
      }
    }
    
  }
  
  void setInitialConditions(state_type &x)
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
        // x[idx2d(i, j)] = initial_value;
        // initial_value *= -1.0; // alternate between -1 and 1 for "random" initial conditions
      }
    }
  }

  double l2residual(const state_type&c){
    state_type dcdt;
    (*this)(c, dcdt, 0);
    double res = 0;
    for (int i = 0; i < nx_*nx_; ++i){
      res += dcdt[i] * dcdt[i];
    }
    return sqrt(res);
  }

private:
  const double D_;     // diffusion coefficient
  const double gamma_; // the other term
  const double b_;
  const double u_;
  const double sig_;
  const double phi_star_;
  const int nx_;       // number of finite difference nodes in each dimension
  const double dx_;       // mesh size

  double laplace_component(double c)
  {
    return D_ * u_ * (c * c * c) - D_ * b_ * c;
  }

  int idx2d_impl(int i, int j)
  {
    return i * nx_ + j;
  }
  
  // regular modulo operator gives negative values without this
  int mod(int a, int b)
  { return (a%b+b)%b; }

  int idx2d(int i, int j)
  {
    // modify the indices to map to a periodic mesh. need two levels for the 4th order operator.
    // i coordinates:
    i = mod(i, nx_);
    j = mod(j, nx_);

    return idx2d_impl(i, j);
  }
};

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



int main()
{
  CHparams chparams;
  
  // *********  Inputs  ***********
  chparams.m        = 1.0;
  chparams.gam      = pow( 0.01 ,2 );
  chparams.b        = 1.0;
  chparams.u        = 1.0;
  chparams.sig      = -10.0;
  chparams.phi_star = 0.2;
  const int nx          = 128;
  const double dx       = 1./nx;
  const int checkpoint  = 50;
  const int maxsteps    = 500;
  // ******************************
  
  CahnHilliard2DRHS rhs(chparams, nx, dx);

  state_type x;
  rhs.setInitialConditions(x);

  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<state_type> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  double time                  = 0.0;
  const double stability_limit = 0.5*dx*dx*dx*dx/chparams.m/chparams.gam; // just an estimate
  double dt_initial            = stability_limit * 0.5;
  double dt_check_residual     = dt_initial * 10.0;
  
  const double res0 = rhs.l2residual(x);

  std::cout << "residual at initial condition: " << res0 << std::endl;
  write_state(x,0,nx);
  
  for (int i = 0; i < maxsteps; ++i){
    integrate_adaptive(controlled_stepper, rhs, x, time, time+dt_check_residual, dt_initial);
    time += dt_check_residual;
    if ( (i+1) % checkpoint == 0) {
      std::cout << "iter: " << i+1 << " , t = " << time << ", relative residual: " << rhs.l2residual(x) / res0 << std::endl;
      write_state(x,i+1,nx);
    }
  }

}
