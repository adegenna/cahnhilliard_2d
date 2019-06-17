#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <omp.h>
#include <boost/numeric/odeint.hpp>
#include "cahnhilliard.h"


int main()
{

  CHparams chparams;
  
  // *********  Inputs  ***********
  chparams.D        = 1.0;
  chparams.gamma    = pow( 0.01 ,2 );
  chparams.b        = 1.0;
  chparams.u        = 1.0;
  chparams.alpha    = 100.0;
  chparams.phi_star = 0.0;
  chparams.sigma    = 10.0;
  chparams.nx       = 128;
  chparams.dx       = 1./chparams.nx;
  chparams.param_type = 0;
  chparams.t0         = 0.0;
  int n_tsteps        = 10;
  double n_dt         = 500.0;
  // ******************************

  double dt_biharm  = (chparams.dx * chparams.dx * chparams.dx * chparams.dx) / chparams.D / chparams.gamma;
  double tf         = n_dt * dt_biharm;
  chparams.dt_check = tf / n_tsteps;

  double dt_diff    = chparams.dx * chparams.dx / chparams.D / chparams.u;
  double dt_lin     = 1.0 / chparams.alpha;

  std::cout << "Biharmonic timescale dt_biharm = " << dt_biharm << std::endl;
  std::cout << "Diffusion timescale dt_diff = " << dt_diff/dt_biharm << " dt_biharm" << std::endl;
  std::cout << "Linear timescale dt_lin = " << dt_lin/dt_biharm << " dt_biharm" << std::endl;

  for (int i=0; i<n_tsteps; i++) {
    chparams.t0 = i * chparams.dt_check;
    chparams.tf = (i+1) * chparams.dt_check;
    std::cout << "t0 = " << chparams.t0/dt_biharm << " dt_biharm , tf = " << chparams.tf/dt_biharm << " dt_biharm" << std::endl;
    run_ch_solver(chparams);
  }
    
}
