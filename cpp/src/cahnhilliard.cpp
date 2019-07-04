#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <omp.h>
#include <boost/numeric/odeint.hpp>
#include <string>
#include "cahnhilliard.h"
#include "cahnhilliard_nonlocal.h"


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
    chpV_.m        = std::vector<double>( info_.nx*info_.nx , chp.m  );
    chpV_.sigma_noise    = chp.sigma_noise;

    if ( info.bc.compare("dirichlet") == 0 ) {
      ch_rhs_ = &compute_ch_nonlocal_stationary_boundaries;
      std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, dirichlet BCs, no thermal dependence" << std::endl;
    }
    else if ( info.bc.compare("neumann") == 0 ) {
      ch_rhs_ = &compute_ch_nonlocal_stationary_boundaries;
      std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, neumann BCs, no thermal dependence" << std::endl;
    }
    else {
      ch_rhs_ = &compute_ch_nonlocal;
      std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, periodic BCs, no thermal dependence" << std::endl;
    }
    
  }

CahnHilliard2DRHS::CahnHilliard2DRHS(CHparamsVector& chp , SimInfo& info)
  : noise_dist_(0.0,1.0) , chpV_(chp) , info_(info)
  {
    if ( info.bc.compare("dirichlet") == 0) {
      ch_rhs_ = &compute_ch_nonlocal_stationary_boundaries;
      std::cout << "Initialized Cahn-Hilliard equation with spatial-field parameters, dirichlet BCs, no thermal dependence" << std::endl;
    }
    else if ( info.bc.compare("neumann") == 0) {
      ch_rhs_ = &compute_ch_nonlocal_stationary_boundaries;
      std::cout << "Initialized Cahn-Hilliard equation with spatial-field parameters, neumann BCs, no thermal dependence" << std::endl;
    }
    else {
      ch_rhs_ = &compute_ch_nonlocal;
      std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, periodic BCs, no thermal dependence" << std::endl;
    }
  }

CahnHilliard2DRHS::~CahnHilliard2DRHS() { };

void CahnHilliard2DRHS::rhs(const std::vector<double> &c, std::vector<double> &dcdt, const double t)
  {
    dcdt.resize(info_.nx*info_.nx);
    (*ch_rhs_)(c, dcdt, t, chpV_, info_);
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
    for (int i = 0; i < info_.nx; ++i) {
      for (int j = 0; j < info_.nx; ++j) {
        x[info_.idx2d(i,j)] = distribution(generator) * 0.005;
      }
    }
    // Set BCs if needed
    if ( info_.bc.compare("dirichlet") == 0) {
      x = apply_dirichlet_bc( x , info_ );
    }
    else if ( info_.bc.compare("neumann") == 0 ) {
      x = apply_neumann_bc( x , info_ );
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


void CahnHilliard2DRHS::write_state(const std::vector<double> &x , const int idx , const int nx )
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
