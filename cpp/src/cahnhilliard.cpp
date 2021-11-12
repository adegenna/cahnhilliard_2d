#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <omp.h>
#include <boost/numeric/odeint.hpp>
#include <string>
#include "cahnhilliard.h"
#include "utils_ch.h"


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
    chpV_.eps_2    = std::vector<double>( info_.nx*info_.ny , chp.eps_2     );
    chpV_.b        = std::vector<double>( info_.nx*info_.ny , chp.b         );
    chpV_.u        = std::vector<double>( info_.nx*info_.ny , chp.u         );
    chpV_.sigma    = std::vector<double>( info_.nx*info_.ny , chp.sigma     );
    chpV_.m        = std::vector<double>( info_.nx*info_.ny , chp.m  );
    chpV_.sigma_noise    = chp.sigma_noise;

    if ( info.bc.compare("dirichlet") == 0 ) {
      ch_rhs_ = &compute_ch_nonlocal_stationary_boundaries;
      std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, dirichlet BCs, no thermal dependence" << std::endl;
    }
    else if ( info.bc.compare("neumann") == 0 ) {
      ch_rhs_ = &compute_ch_nonlocal_neumannBC;
      std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, neumann BCs, no thermal dependence" << std::endl;
    }
    else if ( info.bc.compare("mixed_neumann_bottom_dirichlet") == 0 ) {
      ch_rhs_ = &compute_ch_nonlocal_mixedBC_neumann_with_bottom_dirichlet;
      std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, mixed BC (neumann + bottom dirichlet), no thermal dependence" << std::endl;
    }
    else if ( info.bc.compare("mixed_neumann_top_dirichlet") == 0 ) {
      ch_rhs_ = &compute_ch_nonlocal_mixedBC_neumann_with_top_dirichlet;
      std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, mixed BC (neumann + top dirichlet), no thermal dependence" << std::endl;
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
      ch_rhs_ = &compute_ch_nonlocal_neumannBC;
      std::cout << "Initialized Cahn-Hilliard equation with spatial-field parameters, neumann BCs, no thermal dependence" << std::endl;
    }
    else if ( info.bc.compare("mixed_neumann_bottom_dirichlet") == 0 ) {
      ch_rhs_ = &compute_ch_nonlocal_mixedBC_neumann_with_bottom_dirichlet;
      std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, mixed BC (neumann + bottom dirichlet), no thermal dependence" << std::endl;
    }
    else if ( info.bc.compare("mixed_neumann_top_dirichlet") == 0 ) {
      ch_rhs_ = &compute_ch_nonlocal_mixedBC_neumann_with_top_dirichlet;
      std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, mixed BC (neumann + top dirichlet), no thermal dependence" << std::endl;
    }
    else {
      ch_rhs_ = &compute_ch_nonlocal;
      std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, periodic BCs, no thermal dependence" << std::endl;
    }
  }

CahnHilliard2DRHS::~CahnHilliard2DRHS() { };

void CahnHilliard2DRHS::rhs(const std::vector<double> &c, std::vector<double> &dcdt, const double t)
  {
    dcdt.resize(info_.nx * info_.ny);
    (*ch_rhs_)(c, dcdt, t, chpV_, info_);
  }


void CahnHilliard2DRHS::setInitialConditions(std::vector<double> &x)
  {

    // Set BCs if needed
    if ( info_.bc.compare("dirichlet") == 0) {
      x = apply_dirichlet_bc( x , info_ );
    }
    else if ( info_.bc.compare("mixed_neumann_bottom_dirichlet") == 0 ) {
      x = apply_mixed_bc_neumann_with_bottom_dirichlet( x , info_ );
    }
    else if ( info_.bc.compare("mixed_neumann_top_dirichlet") == 0 ) {
      x = apply_mixed_bc_neumann_with_top_dirichlet( x , info_ );
    }
    else if ( info_.bc.compare("neumann") == 0 ) {
      x = apply_neumann_bc( x , info_ );
    }
    
  }

void CahnHilliard2DRHS::write_state(const std::vector<double> &x , const int idx , const int nx , const int ny , std::string& outdir)
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
