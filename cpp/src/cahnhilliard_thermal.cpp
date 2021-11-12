#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <omp.h>
#include <boost/numeric/odeint.hpp>
#include "cahnhilliard_thermal.h"
#include "utils_ch.h"

  /*
  Cahn-Hilliard:
  
  dc/dt = laplacian( u*c^3 - b*c ) - eps_2*biharm(c) - sigma*(c - m) + sigma_noise * N(0,1^2)
  
  expanding out RHS into individual differentials:
  D*laplacian( u*c^3 - b*c) - D*eps_2*biharm(c)
  assuming constant eps_2.

  need a d^4 and a d^2 operator.
  */

CahnHilliard2DRHS_thermal::CahnHilliard2DRHS_thermal(CHparamsScalar& chp , SimInfo& info)
  : noise_dist_(0.0,1.0) , info_(info)
  {    
    chpV_.eps_2    = std::vector<double>( info_.nx*info_.ny , chp.eps_2     );
    chpV_.b        = std::vector<double>( info_.nx*info_.ny , chp.b         );
    chpV_.u        = std::vector<double>( info_.nx*info_.ny , chp.u         );
    chpV_.sigma    = std::vector<double>( info_.nx*info_.ny , chp.sigma     );
    chpV_.m        = std::vector<double>( info_.nx*info_.ny , chp.m  );
    chpV_.DT       = std::vector<double>( info_.nx*info_.ny , chp.DT  );
    chpV_.f_T      = std::vector<double>( info_.nx*info_.ny , chp.f_T  );
    chpV_.sigma_noise    = chp.sigma_noise;

    if ( info.bc.compare("dirichlet") == 0) {
      ch_rhs_ = &compute_ch_nonlocal_stationary_boundaries;
      std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, dirichlet BCs, thermal coefficient dependence, thermal diffusion" << std::endl;
    }
    else if ( info.bc.compare("neumann") == 0) {
      ch_rhs_ = &compute_ch_nonlocal_neumannBC;
      std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, neumann BCs, thermal coefficient dependence, thermal diffusion" << std::endl;
    }
    else {
      ch_rhs_ = &compute_ch_nonlocal;
      std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, periodic BCs, thermal coefficient dependence, thermal diffusion" << std::endl;
    }
    
  }

CahnHilliard2DRHS_thermal::CahnHilliard2DRHS_thermal(CHparamsVector& chp , SimInfo& info)
  : noise_dist_(0.0,1.0) , chpV_(chp) , info_(info)
  {

    if ( info.bc.compare("dirichlet") == 0) {
      ch_rhs_ = &compute_ch_nonlocal_stationary_boundaries;
      std::cout << "Initialized Cahn-Hilliard equation: spatial-field parameters, dirichlet BCs, thermal coefficient dependence, thermal diffusion" << std::endl;
    }
    if ( info.bc.compare("neumann") == 0) {
      ch_rhs_ = &compute_ch_nonlocal_neumannBC;
      std::cout << "Initialized Cahn-Hilliard equation: spatial-field parameters, neumann BCs, thermal coefficient dependence, thermal diffusion" << std::endl;
    }
    else {
      ch_rhs_ = &compute_ch_nonlocal;
      std::cout << "Initialized Cahn-Hilliard equation: spatial-field parameters, periodic BCs, thermal coefficient dependence, thermal diffusion" << std::endl;
    }
    
  }

CahnHilliard2DRHS_thermal::~CahnHilliard2DRHS_thermal() { };

void CahnHilliard2DRHS_thermal::rhs(const std::vector<double> &ct, std::vector<double> &dcTdt, const double t)
  {
    dcTdt.resize(2 * info_.nx * info_.ny);
    std::vector<double> c = std::vector<double>( ct.begin() , ct.begin() + info_.nx*info_.ny );
    std::vector<double> T = std::vector<double>( ct.begin() + info_.nx*info_.ny , ct.end() );

    // enforce thermal BC: dT/dnormal = 0
    # pragma omp parallel for
    for (int i = 0; i < info_.nx; ++i) {
      T[info_.idx2d(0,i)]          = T[info_.idx2d(1,i)];
      T[info_.idx2d(info_.ny-1,i)] = T[info_.idx2d(info_.ny-2,i)];
    }

    # pragma omp parallel for
    for (int i = 0; i < info_.ny; ++i) {
      T[info_.idx2d(i,0)]          = T[info_.idx2d(i,1)];
      T[info_.idx2d(i,info_.nx-1)] = T[info_.idx2d(i,info_.nx-2)];
    }

    // evaluate CH parameter dependencies on temperature
    //chpV_ = compute_chparams_using_temperature( chpV_ , info_ , T );
    chpV_ = compute_eps2_and_sigma_from_polymer_params( chpV_ , info_ , T );

    // evaluate deterministic nonlocal dynamics
    (*ch_rhs_)(c, dcTdt, t, chpV_, info_);
    
    // evaluate thermal diffusion
    # pragma omp parallel for
    for (int i = 0; i < info_.ny; ++i) {
      for (int j = 0; j < info_.nx; ++j) {
        
        const double T_i   = T[info_.idx2d(i, j)];
        const double T_im1 = T[info_.idx2d(i - 1, j)];
        const double T_ip1 = T[info_.idx2d(i + 1, j)];
        const double T_jm1 = T[info_.idx2d(i, j - 1)];
        const double T_jp1 = T[info_.idx2d(i, j + 1)];

        double dxx = 1.0 / (info_.dx * info_.dx) * (T_jm1 + T_jp1 - 2.0 * T_i);
        double dyy = 1.0 / (info_.dy * info_.dy) * (T_im1 + T_ip1 - 2.0 * T_i);
        
        dcTdt[info_.idx2d(i, j) + info_.nx*info_.ny]  = chpV_.DT[info_.idx2d(i, j)] * (dxx + dyy) + chpV_.f_T[info_.idx2d(i, j)];
        
      }
    }

    // enforce thermal BC: dT/dnormal = 0
    # pragma omp parallel for
    for (int i = 0; i < info_.nx; ++i) {
      dcTdt[info_.idx2d(0,i) + info_.nx*info_.ny]          = 0;
      dcTdt[info_.idx2d(info_.ny-1,i) + info_.nx*info_.ny] = 0;
    }

    # pragma omp parallel for
    for (int i = 0; i < info_.ny; ++i) {
      dcTdt[info_.idx2d(i,0) + info_.nx*info_.ny]          = 0;
      dcTdt[info_.idx2d(i,info_.nx-1) + info_.nx*info_.ny] = 0;
    }
    
  }


void CahnHilliard2DRHS_thermal::setInitialConditions(std::vector<double> &x)
  {
    
    // Set BCs if needed
    if ( info_.bc.compare("dirichlet") == 0) {
      x = apply_dirichlet_bc( x , info_ );
    }
    else if ( info_.bc.compare("neumann") == 0 ) {
      x = apply_neumann_bc( x , info_ );
    }

  }


void CahnHilliard2DRHS_thermal::write_state(const std::vector<double> &x , const int idx , const int nx , const int ny , std::string& outdir)
{
  if ( outdir.back() != '/' )
    outdir += '/';
  std::ofstream outC;
  std::ofstream outT;
  outC.open( outdir + "C_" + std::to_string(idx) + ".out" );
  outT.open( outdir + "T_" + std::to_string(idx) + ".out" );
  outC.precision(16);
  outT.precision(16);
  
  for (int i = 0; i < ny; ++i){
    for (int j = 0; j < nx; ++j){
      outC << x[i * ny + j] << " ";
      outT << x[i * ny + j + nx*ny] << " ";
    }
  }

  outC.close();
  outT.close();
};
