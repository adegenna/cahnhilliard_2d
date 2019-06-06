#include <iostream>
#include <vector>
#include <random>

#include <boost/numeric/odeint.hpp>

#include "matplotlibcpp.h"

typedef std::vector<double> state_type;

class CahnHilliard2DRHS
{
public:
  CahnHilliard2DRHS(double D, double g, int nx, double dx)
      : D_(D), gamma_(g), nx_(nx), dx_(dx)
  {
    std::cout << "Initialized Cahn-Hilliard equation with D_ " << D_ 
      << " gamma_ " << gamma_ << " dx_ " << nx_ << " dx_ " << dx_ << std::endl;
  }

  /*
  Cahn-Hilliard according to https://en.wikipedia.org/wiki/Cahn%E2%80%93Hilliard_equation is
  dc/dt = D laplacian(c^3 - c - gamma*laplacian(c))
  expanding out RHS into individual differentials:
  D laplacian(c^3 - c) - D gamma laplacian(laplacian(c))
  assuming constant gamma.

  need a d^4 and a d^2 operator.
  */
  void operator()(const state_type &c, state_type &dcdt, const double t)
  {
    dcdt.resize(nx_*nx_);

    // evaluate the second order term, 5 point central stencil
    for (int i = 0; i < nx_; ++i)
    {
      for (int j = 0; j < nx_; ++j)
      {
        const double c_i = laplace_component(c[idx2d(i, j)]);
        const double c_im1 = laplace_component(c[idx2d(i - 1, j)]);
        const double c_ip1 = laplace_component(c[idx2d(i + 1, j)]);
        const double c_jm1 = laplace_component(c[idx2d(i, j - 1)]);
        const double c_jp1 = laplace_component(c[idx2d(i, j + 1)]);
        dcdt[idx2d(i, j)] = (D_ / (dx_ * dx_)) * (c_im1 + c_ip1 + c_jm1 + c_jp1 - 4.0 * c_i);
      }
    }

    // evaluate the 4th order term, 9 point central stencil
    for (int i = 0; i < nx_; ++i){
      for (int j = 0; j < nx_; ++j){
        const double c_i = c[idx2d(i, j)];
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
        x[idx2d(i,j)] = distribution(generator);
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
  const int nx_;       // number of finite difference nodes in each dimension
  const double dx_;       // mesh size

  double laplace_component(double c)
  {
    return c * c * c - c;
    // return c;
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

int main()
{
  const double D = 1.0;
  const double gam = 0.5;
  const int nx = 128;
  const double dx = 0.1;

  CahnHilliard2DRHS rhs(D, gam, nx, dx);

  state_type x;
  rhs.setInitialConditions(x);

  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<state_type> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  double time = 0.0;
  const double stability_limit = 0.5*dx*dx*dx*dx/D; // just an estimate
  double dt_initial = stability_limit * 0.5;
  double dt_check_residual = dt_initial * 10.0;
  const int maxsteps = 50;

  const double res0 = rhs.l2residual(x);

  std::cout << "residual at initial condition: " << res0 << std::endl;

  for (int i = 0; i < maxsteps; ++i){
    integrate_adaptive(controlled_stepper, rhs, x, time, time+dt_check_residual, dt_initial);
    std::cout << "iter: " << i << " relative residual: " << rhs.l2residual(x) / res0 << std::endl;
    time += dt_check_residual;
  }


  std::vector<std::vector<double>> xx, yy, zz;
  for (int i = 0; i < nx; ++i){
    std::vector<double> x_row, y_row, z_row;
    for (int j = 0; j < nx; ++j){
      x_row.push_back(i*dx);
      y_row.push_back(j*dx);
      z_row.push_back(x[i * nx + j]);
    }
    xx.push_back(x_row);
    yy.push_back(y_row);
    zz.push_back(z_row);
  }

  namespace plt = matplotlibcpp;
  plt::plot_surface(xx, yy, zz);
  plt::show();
}
