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
  chparams.m        = 1.0;
  chparams.gam      = pow( 0.01 ,2 );
  chparams.b        = 1.0;
  chparams.u        = 1.0;
  chparams.alpha    = 10.0;
  chparams.phi_star = 0.0;
  chparams.sigma    = 0.0;
  const int nx          = 128;
  const double dx       = 1./nx;
  const int checkpoint  = 20;
  const int maxsteps    = 200;
  // ******************************

  run_ch_solver(chparams, nx, dx, checkpoint, maxsteps);
  
}
