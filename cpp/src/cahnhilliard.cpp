#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <omp.h>
#include <boost/numeric/odeint.hpp>
#include "cahnhilliard.h"
//#include "stochastic_euler.hpp"


  /*
  Cahn-Hilliard:
  
  dc/dt = D*laplacian( u*c^3 - b*c) - D*gamma*biharm(c)
  
  expanding out RHS into individual differentials:
  D*laplacian( u*c^3 - b*c) - D*gamma*biharm(c)
  assuming constant gamma.

  need a d^4 and a d^2 operator.
  */


CahnHilliard2DRHS::CahnHilliard2DRHS(CHparamsScalar& chp , SimInfo& info)
  : noise_dist_(0.0,1.0) , chpS_(chp) , info_(info)
  {
    std::cout << "Initialized Cahn-Hilliard equation with scalar parameters" << std::endl;
  }

CahnHilliard2DRHS::CahnHilliard2DRHS(CHparamsVector& chp , SimInfo& info)
  : noise_dist_(0.0,1.0) , chpV_(chp) , info_(info)
  {
    std::cout << "Initialized Cahn-Hilliard equation with spatial-field parameters" << std::endl;
  }

CahnHilliard2DRHS::~CahnHilliard2DRHS() { };

void CahnHilliard2DRHS_Scalar::get_ij_values(int i , int j , CHparamsScalar& ch_ij) {
  ch_ij.D        = chpS_.D;
  ch_ij.gamma    = chpS_.gamma;
  ch_ij.b        = chpS_.b;
  ch_ij.u        = chpS_.u;
  ch_ij.alpha    = chpS_.alpha;
  ch_ij.phi_star = chpS_.phi_star;
  ch_ij.sigma    = chpS_.sigma;
  
};


void CahnHilliard2DRHS_Vector::get_ij_values(int i , int j , CHparamsScalar& ch_ij) {
  ch_ij.D        = chpV_.D[idx2d(i, j)];
  ch_ij.gamma    = chpV_.gamma[idx2d(i, j)];
  ch_ij.b        = chpV_.b[idx2d(i, j)];
  ch_ij.u        = chpV_.u[idx2d(i, j)];
  ch_ij.alpha    = chpV_.alpha[idx2d(i, j)];
  ch_ij.phi_star = chpV_.phi_star[idx2d(i, j)];
  ch_ij.sigma    = chpV_.sigma;
  
};


void CahnHilliard2DRHS::rhs(const std::vector<double> &c, std::vector<double> &dcdt, const double t)
  {
    dcdt.resize(info_.nx*info_.nx);

    // evaluate the second order term, 5 point central stencil
    # pragma omp parallel for
    for (int i = 0; i < info_.nx; ++i)
    {
      for (int j = 0; j < info_.nx; ++j)
      {
        get_ij_values(i,j,ch_ij_);
        
        const double c_i   = laplace_component(c[idx2d(i, j)]     , ch_ij_.D , ch_ij_.u , ch_ij_.b);
        const double c_im1 = laplace_component(c[idx2d(i - 1, j)] , ch_ij_.D , ch_ij_.u , ch_ij_.b);
        const double c_ip1 = laplace_component(c[idx2d(i + 1, j)] , ch_ij_.D , ch_ij_.u , ch_ij_.b);
        const double c_jm1 = laplace_component(c[idx2d(i, j - 1)] , ch_ij_.D , ch_ij_.u , ch_ij_.b);
        const double c_jp1 = laplace_component(c[idx2d(i, j + 1)] , ch_ij_.D , ch_ij_.u , ch_ij_.b);
        
        dcdt[idx2d(i, j)]  = (ch_ij_.D / (info_.dx * info_.dx)) * (c_im1 + c_ip1 + c_jm1 + c_jp1 - 4.0 * c_i);
      }
    }

    // evaluate the 4th order term, 9 point central stencil
    # pragma omp parallel for
    for (int i = 0; i < info_.nx; ++i){
      for (int j = 0; j < info_.nx; ++j){
        get_ij_values(i,j,ch_ij_);
        
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
        dcdt[idx2d(i,j)] -= ch_ij_.D * ch_ij_.gamma /(info_.dx*info_.dx*info_.dx*info_.dx) * 
          (c_ip2 - 4.0*c_ip1 + 6.0*c_i - 4.0*c_im1 + c_im2);

        // y-direction u_yyyy
        dcdt[idx2d(i,j)] -= ch_ij_.D * ch_ij_.gamma /(info_.dx*info_.dx*info_.dx*info_.dx) * 
          (c_jp2 - 4.0*c_jp1 + 6.0*c_i - 4.0*c_jm1 + c_jm2);

        // mixed term 2*u_xxyy
        dcdt[idx2d(i,j)] -= ch_ij_.D * ch_ij_.gamma /(info_.dx*info_.dx*info_.dx*info_.dx) * 
          2 * (4*c_i - 2*(c_im1 + c_ip1 + c_jm1 + c_jp1) + c_ul + c_ur + c_bl + c_br );
      }
    }

    // evaluate linear term
    # pragma omp parallel for
    for (int i = 0; i < info_.nx; ++i){
      for (int j = 0; j < info_.nx; ++j){
        get_ij_values(i,j,ch_ij_);
        
        const double c_i   = c[idx2d(i, j)];
        
        dcdt[idx2d(i,j)]  -= ch_ij_.alpha * ( c_i - ch_ij_.phi_star );
      }
    }
    
  }


void CahnHilliard2DRHS::operator()(const std::vector<double> &c, std::vector<double> &dcdt, const double t)
{
  rhs(c,dcdt,t);
}

void CahnHilliard2DRHS::setInitialConditions(std::vector<double> &x)
  {
    x.resize(info_.nx * info_.nx);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-1.0,1.0);

    // double initial_value = -1.0;
    for (int i = 0; i < info_.nx; ++i)
    {
      for (int j = 0; j < info_.nx; ++j)
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
    for (int i = 0; i < info_.nx*info_.nx; ++i){
      res += dcdt[i] * dcdt[i];
    }
    return sqrt(res);
  }

double CahnHilliard2DRHS::laplace_component(double c , double D , double u , double b)
  {
    return D * u * (c * c * c) - D * b * c;
  }

int CahnHilliard2DRHS::idx2d_impl(int i, int j)
  {
    return i * info_.nx + j;
  }
  
  // regular modulo operator gives negative values without this
int CahnHilliard2DRHS::mod(int a, int b)
  { return (a%b+b)%b; }

int CahnHilliard2DRHS::idx2d(int i, int j)
  {
    // modify the indices to map to a periodic mesh. need two levels for the 4th order operator.
    // i coordinates:
    i = mod(i, info_.nx);
    j = mod(j, info_.nx);

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

// void run_ch_solver_checkpointing(CHparams& chparams , SimInfo& info)
// {
//   CahnHilliard2DRHS rhs(chparams);

//   std::vector<double> x;
//   rhs.setInitialConditions(x);

//   // define adaptive stepper
//   typedef boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>> error_stepper_type;

//   // define runge kutta
//   typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

//   controlled_stepper_type controlled_stepper;

//   double time                  = 0.0;
//   const double stability_limit = 0.5*info.dx*info.dx*info.dx*info.dx/chparams.D/chparams.gamma; // just an estimate
//   double dt_initial            = stability_limit * 0.5;
  
//   const double res0            = rhs.l2residual(x);

//   std::cout << "residual at initial condition: " << res0 << std::endl;
//   write_state(x,0,chparams.nx);

//   double t_steps = int( chparams.tf / chparams.dt_check );
//   for (int i = 0; i < t_steps; ++i){
//     integrate_adaptive(controlled_stepper, rhs, x, time, time + chparams.dt_check, dt_initial);
//     time += chparams.dt_check;
//     std::cout << "iter: " << i+1 << " , t = " << time << ", relative residual: " << rhs.l2residual(x) / res0 << std::endl;
//     write_state(x,i+1,chparams.nx);
//   }

// };

