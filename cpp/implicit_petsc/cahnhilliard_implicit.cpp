#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <omp.h>
#include <boost/numeric/odeint.hpp>
#include <string>
#include "cahnhilliard_implicit.h"
#include "utils_ch.h"


  /*
  Cahn-Hilliard:
  
  dc/dt = laplacian( u*c^3 - b*c ) - eps_2*biharm(c) - sigma*(c - m) + sigma_noise * N(0,1^2)
  
  expanding out RHS into individual differentials:
  D*laplacian( u*c^3 - b*c) - D*eps_2*biharm(c)
  assuming constant eps_2.

  need a d^4 and a d^2 operator.
  */

CahnHilliardImplicit2DRHS::CahnHilliardImplicit2DRHS(CHparamsScalar& chp , SimInfo& info , DM& da , TS& ts , Vec& F)
  : noise_dist_(0.0,1.0) , info_(info) , da_(da) , ts_(ts) , F_(F)
  {    
    chpV_.eps_2    = std::vector<double>( info_.nx*info_.ny , chp.eps_2     );
    chpV_.b        = std::vector<double>( info_.nx*info_.ny , chp.b         );
    chpV_.u        = std::vector<double>( info_.nx*info_.ny , chp.u         );
    chpV_.sigma    = std::vector<double>( info_.nx*info_.ny , chp.sigma     );
    chpV_.m        = std::vector<double>( info_.nx*info_.ny , chp.m  );
    chpV_.sigma_noise    = chp.sigma_noise;
    
    if ( info.bc.compare("dirichlet") == 0 ) {
      ch_rhs_ = &compute_ch_nonlocal_implicit_dirichlet;
      std::cout << "Initialized Cahn-Hilliard equation: implicit, scalar parameters, dirichlet BCs, no thermal dependence" << std::endl;
    }
    
  }

CahnHilliardImplicit2DRHS::CahnHilliardImplicit2DRHS(CHparamsVector& chp , SimInfo& info)
  : noise_dist_(0.0,1.0) , chpV_(chp) , info_(info)
  {
    if ( info.bc.compare("dirichlet") == 0) {
      ch_rhs_ = &compute_ch_nonlocal_implicit_dirichlet;
      std::cout << "Initialized Cahn-Hilliard equation: implicit, spatial-field parameters, dirichlet BCs, no thermal dependence" << std::endl;
    }
  }

CahnHilliardImplicit2DRHS::~CahnHilliardImplicit2DRHS() { };

void CahnHilliardImplicit2DRHS::rhs(const std::vector<double> &c, std::vector<double> &dcdt, const double t)
  {
    dcdt.resize(info_.nx * info_.ny);
    (*ch_rhs_)(c, dcdt, t, chpV_, info_, ts_ , F_ , da_);
  }


void CahnHilliardImplicit2DRHS::setInitialConditions(std::vector<double> &x)
  {
    x.resize(info_.nx * info_.ny);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-1.0,1.0);

    // double initial_value = -1.0;
    for (int i = 0; i < info_.ny; ++i) {
      for (int j = 0; j < info_.nx; ++j) {
        x[info_.idx2d(i,j)] = distribution(generator) * 0.005;
      }
    }
    
  }

void CahnHilliardImplicit2DRHS::write_state(const std::vector<double> &x , const int idx , const int nx , const int ny , std::string& outdir)
{
  if ( outdir.back() != '/' )
    outdir += '/';
  std::ofstream out;
  out.open( outdir + "C_" + std::to_string(idx) + ".out" ); 
  out.precision(16);
  
  for (int i = 0; i < nx; ++i){
    for (int j = 0; j < ny; ++j){
      out << x[i * ny + j] << " ";
    }
  }

  out.close();
};
