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
    chpV_.eps_2    = std::vector<double>( info_.nx*info_.nx , chp.eps_2     );
    chpV_.b        = std::vector<double>( info_.nx*info_.nx , chp.b         );
    chpV_.u        = std::vector<double>( info_.nx*info_.nx , chp.u         );
    chpV_.sigma    = std::vector<double>( info_.nx*info_.nx , chp.sigma     );
    chpV_.m        = std::vector<double>( info_.nx*info_.nx , chp.m  );
    chpV_.DT       = std::vector<double>( info_.nx*info_.nx , chp.DT  );
    chpV_.f_T      = std::vector<double>( info_.nx*info_.nx , chp.f_T  );
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
    dcTdt.resize(2*info_.nx*info_.nx);
    std::vector<double> c = std::vector<double>( ct.begin() , ct.begin() + info_.nx*info_.nx );
    std::vector<double> T = std::vector<double>( ct.begin() + info_.nx*info_.nx , ct.end() );

    // enforce thermal BC: dT/dnormal = 0
    # pragma omp parallel for
    for (int i = 0; i < info_.nx; ++i) {
      T[info_.idx2d(i,0)]          = T[info_.idx2d(i,1)];
      T[info_.idx2d(i,info_.nx-1)] = T[info_.idx2d(i,info_.nx-2)];
      T[info_.idx2d(0,i)]          = T[info_.idx2d(1,i)];
      T[info_.idx2d(info_.nx-1,i)] = T[info_.idx2d(info_.nx-2,i)];
    }

    // evaluate CH parameter dependencies on temperature
    chpV_ = compute_chparams_using_temperature( chpV_ , info_ , T );

    // evaluate deterministic nonlocal dynamics
    compute_ch_nonlocal(c, dcTdt, t, chpV_, info_);
    
    // evaluate thermal diffusion
    # pragma omp parallel for
    for (int i = 0; i < info_.nx; ++i) {
      for (int j = 0; j < info_.nx; ++j) {
        
        const double T_i   = T[info_.idx2d(i, j)];
        const double T_im1 = T[info_.idx2d(i - 1, j)];
        const double T_ip1 = T[info_.idx2d(i + 1, j)];
        const double T_jm1 = T[info_.idx2d(i, j - 1)];
        const double T_jp1 = T[info_.idx2d(i, j + 1)];

        dcTdt[info_.idx2d(i, j) + info_.nx*info_.nx]  = (chpV_.DT[info_.idx2d(i, j)] / (info_.dx * info_.dx)) * (T_im1 + T_ip1 + T_jm1 + T_jp1 - 4.0 * T_i) + chpV_.f_T[info_.idx2d(i, j)];
      }
    }

    // enforce thermal BC: dT/dnormal = 0
    # pragma omp parallel for
    for (int i = 0; i < info_.nx; ++i) {
      dcTdt[info_.idx2d(i,0) + info_.nx*info_.nx]          = 0;
      dcTdt[info_.idx2d(i,info_.nx-1) + info_.nx*info_.nx] = 0;
      dcTdt[info_.idx2d(0,i) + info_.nx*info_.nx]          = 0;
      dcTdt[info_.idx2d(info_.nx-1,i) + info_.nx*info_.nx] = 0;
    }
    
  }


void CahnHilliard2DRHS_thermal::operator()(const std::vector<double> &c, std::vector<double> &dcdt, const double t)
{
  rhs(c,dcdt,t);
}

void CahnHilliard2DRHS_thermal::setInitialConditions(std::vector<double> &x)
  {
    x.resize(2 * info_.nx * info_.nx);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-1.0,1.0);

    for (int i = 0; i < info_.nx; ++i) {
      for (int j = 0; j < info_.nx; ++j) {
        x[info_.idx2d(i,j)]                     = distribution(generator) * 0.005;
	x[info_.idx2d(i,j) + info_.nx*info_.nx] = chpV_.T_min;
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

double CahnHilliard2DRHS_thermal::l2residual(const std::vector<double>&cT)
  {
    std::vector<double> dcTdt;
    (*this)(cT, dcTdt, 0);
    double res = 0;
    for (int i = 0; i < 2*info_.nx*info_.nx; ++i){
      res += dcTdt[i] * dcTdt[i];
    }
    return sqrt(res);
  }

void CahnHilliard2DRHS_thermal::write_state(const std::vector<double> &x , const int idx , const int nx )
{
  std::cout << "printing C and T ..." << std::endl;
  std::ofstream outC;
  std::ofstream outT;
  outC.open( "C_" + std::to_string(idx) + ".out" );
  outT.open( "T_" + std::to_string(idx) + ".out" );
  outC.precision(16);
  outT.precision(16);
  
  for (int i = 0; i < nx; ++i){
    for (int j = 0; j < nx; ++j){
      outC << x[i * nx + j] << " ";
      outT << x[i * nx + j + nx*nx] << " ";
    }
  }

  outC.close();
  outT.close();
};
