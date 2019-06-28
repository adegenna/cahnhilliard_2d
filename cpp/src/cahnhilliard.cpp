#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <omp.h>
#include <boost/numeric/odeint.hpp>
#include "cahnhilliard.h"


  /*
  Cahn-Hilliard:
  
  dc/dt = laplacian( u*c^3 - b*c ) - eps_2*biharm(c) - sigma*(c - m) + sigma_noise * N(0,1^2)
  
  expanding out RHS into individual differentials:
  D*laplacian( u*c^3 - b*c) - D*eps_2*biharm(c)
  assuming constant eps_2.

  need a d^4 and a d^2 operator.
  */

CahnHilliard2DRHS::CahnHilliard2DRHS(CHparamsScalar& chp , SimInfo& info)
  : noise_dist_(0.0,1.0) , info_(info)
  {    
    chpV_.eps_2    = std::vector<double>( info_.nx*info_.nx , chp.eps_2     );
    chpV_.b        = std::vector<double>( info_.nx*info_.nx , chp.b         );
    chpV_.u        = std::vector<double>( info_.nx*info_.nx , chp.u         );
    chpV_.sigma    = std::vector<double>( info_.nx*info_.nx , chp.sigma     );
    chpV_.m = std::vector<double>( info_.nx*info_.nx , chp.m  );
    chpV_.sigma_noise    = chp.sigma_noise;
    
    std::cout << "Initialized Cahn-Hilliard equation with scalar parameters" << std::endl;
  }

CahnHilliard2DRHS::CahnHilliard2DRHS(CHparamsVector& chp , SimInfo& info)
  : noise_dist_(0.0,1.0) , chpV_(chp) , info_(info)
  {
    std::cout << "Initialized Cahn-Hilliard equation with spatial-field parameters" << std::endl;
  }

CahnHilliard2DRHS::~CahnHilliard2DRHS() { };

void CahnHilliard2DRHS::rhs(const std::vector<double> &c, std::vector<double> &dcdt, const double t)
  {
    dcdt.resize(info_.nx*info_.nx);

    // evaluate the second order term, 5 point central stencil
    # pragma omp parallel for
    for (int i = 0; i < info_.nx; ++i) {
      for (int j = 0; j < info_.nx; ++j) {
        
        const double c_i   = laplace_component( idx2d(i, j)      , c , chpV_.u , chpV_.b );
        const double c_im1 = laplace_component( idx2d(i - 1, j)  , c , chpV_.u , chpV_.b );
        const double c_ip1 = laplace_component( idx2d(i + 1, j)  , c , chpV_.u , chpV_.b );
        const double c_jm1 = laplace_component( idx2d(i, j - 1)  , c , chpV_.u , chpV_.b );
        const double c_jp1 = laplace_component( idx2d(i, j + 1)  , c , chpV_.u , chpV_.b );
        
        dcdt[idx2d(i, j)]  = (1.0 / (info_.dx * info_.dx)) * (c_im1 + c_ip1 + c_jm1 + c_jp1 - 4.0 * c_i);
      }
    }

    // evaluate the 4th order term, 9 point central stencil
    # pragma omp parallel for
    for (int i = 0; i < info_.nx; ++i) {
      for (int j = 0; j < info_.nx; ++j) {
        
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
        dcdt[idx2d(i,j)] -= chpV_.eps_2[idx2d(i,j)] /(info_.dx*info_.dx*info_.dx*info_.dx) * 
          (c_ip2 - 4.0*c_ip1 + 6.0*c_i - 4.0*c_im1 + c_im2);

        // y-direction u_yyyy
        dcdt[idx2d(i,j)] -= chpV_.eps_2[idx2d(i,j)] /(info_.dx*info_.dx*info_.dx*info_.dx) * 
          (c_jp2 - 4.0*c_jp1 + 6.0*c_i - 4.0*c_jm1 + c_jm2);

        // mixed term 2*u_xxyy
        dcdt[idx2d(i,j)] -= chpV_.eps_2[idx2d(i,j)] /(info_.dx*info_.dx*info_.dx*info_.dx) * 
          2 * (4*c_i - 2*(c_im1 + c_ip1 + c_jm1 + c_jp1) + c_ul + c_ur + c_bl + c_br );
      }
    }

    // evaluate linear term
    # pragma omp parallel for
    for (int i = 0; i < info_.nx; ++i){
      for (int j = 0; j < info_.nx; ++j){
        
        const double c_i   = c[idx2d(i, j)];
        
        dcdt[idx2d(i,j)]  -= chpV_.sigma[idx2d(i,j)] * ( c_i - chpV_.m[idx2d(i,j)] );
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

double CahnHilliard2DRHS::laplace_component(int i ,
                           const std::vector<double>& c ,
			   const std::vector<double>& u ,
                           const std::vector<double>& b )
  {
    return u[i] * (c[i] * c[i] * c[i]) - b[i] * c[i];
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
