#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <omp.h>
#include <boost/numeric/odeint.hpp>
#include "cahnhilliard.h"


int main()
{

  CHparamsScalar chparams;
  SimInfo info;
  
  // *********  Inputs  ***********
  chparams.D        = 1.0;
  chparams.gamma    = pow( 0.01 ,2 );
  chparams.b        = 1.0;
  chparams.u        = 1.0;
  chparams.alpha    = 100.0;
  chparams.phi_star = 0.0;
  chparams.sigma    = 10.0;

  info.nx             = 128;
  info.dx             = 1./info.nx;
  info.t0             = 0.0;
  int n_tsteps        = 10;
  double n_dt         = 500.0;
  // ******************************

  double dt_biharm  = (info.dx * info.dx * info.dx * info.dx) / chparams.D / chparams.gamma;
  double tf         = n_dt * dt_biharm;
  info.dt_check     = tf / n_tsteps;

  double dt_diff    = info.dx * info.dx / chparams.D / chparams.u;
  double dt_lin     = 1.0 / chparams.alpha;

  std::cout << "Biharmonic timescale dt_biharm = " << dt_biharm << std::endl;
  std::cout << "Diffusion timescale dt_diff = " << dt_diff/dt_biharm << " dt_biharm" << std::endl;
  std::cout << "Linear timescale dt_lin = " << dt_lin/dt_biharm << " dt_biharm" << std::endl;

  for (int i=0; i<n_tsteps; i++) {
    info.t0 = i * info.dt_check;
    info.tf = (i+1) * info.dt_check;
    std::cout << "t0 = " << info.t0/dt_biharm << " dt_biharm , tf = " << info.tf/dt_biharm << " dt_biharm" << std::endl;
    run_ch_solver(chparams , info);
  }
    
}
