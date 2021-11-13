#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <omp.h>
#include "run_ch_solver.h"

int main()
{

  // ***************************************************************
  // Example driver program for 2D modified Cahn-Hilliard
  // No thermal dependence
  // No polymer property dependence
  // No thermal dynamics
  // Scalar coefficients
  // ***************************************************************

  CHparamsScalar chparams;
  
  int nx=128; int ny=128;
  std::vector<double> x0( nx*ny );
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(-1.0,1.0);
  for (int i = 0; i < ny; ++i) {
    for (int j = 0; j < nx; ++j) {
      x0[ i*nx + j ] = distribution(generator) * 0.005;
    }
  }

  SimInfo info( 128 , 128 , x0 );

  // *********  Inputs  ***********
  double eps_2          = pow( 0.01 ,2 );
  chparams.eps_2        = eps_2;
  chparams.b            = eps_2 / info.dx / info.dx;
  chparams.u            = eps_2 / info.dx / info.dx;
  chparams.sigma        = eps_2 / info.dx / info.dx / info.dx / info.dx / 200.0;
  chparams.m            = 0.0;
  chparams.sigma_noise  = 0.0;
  
  int n_tsteps        = 25;
  double n_dt         = 300.0;
  // ******************************

  // Compute linear timescales
  double dt_biharm  = (info.dx * info.dx * info.dx * info.dx) / chparams.eps_2;
  double dt_diff    = info.dx * info.dx / chparams.u;
  double dt_lin     = 1.0 / chparams.sigma;

  // Setup temporal checkpointing
  double tf         = n_dt * dt_biharm;
  double dt_check   = tf / n_tsteps;


  std::cout << "Biharmonic timescale dt_biharm = " << dt_biharm << std::endl;
  std::cout << "Diffusion timescale dt_diff = " << dt_diff/dt_biharm << " dt_biharm" << std::endl;
  std::cout << "Linear timescale dt_lin = " << dt_lin/dt_biharm << " dt_biharm" << std::endl;

  // Run solver
  for (int i=0; i<n_tsteps; i++) {
    info.t0 = i * dt_check;
    info.tf = (i+1) * dt_check;
    std::cout << "t0 = " << info.t0/dt_biharm << " dt_biharm , tf = " << info.tf/dt_biharm << " dt_biharm" << std::endl;
    run_ch_solver(chparams , info);
  }
    
}
