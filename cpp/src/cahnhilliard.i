 %module cahnhilliard
 %{
 /* Put header files here or function declarations like below */
 #include <iostream>
 #include <vector>
 #include <random>
 #include <fstream>
 #include <omp.h>
 #include <boost/numeric/odeint.hpp>
 #include "right_hand_side.h"
 #include "cahnhilliard.h"
 #include "cahnhilliard_thermal.h"
 #include "cahnhilliard_thermal_nodiffusion.h"
 #include "chparams.h"
 #include "run_ch_solver.hpp"
 #include "utils_ch.h"
 %}

 %include <std_vector.i>
 %include typemaps.i
 %include std_string.i
 namespace std {
  %template(DoubleVector) vector<double>;
 };

 %include "right_hand_side.h"
 %include "cahnhilliard.h"
 %include "cahnhilliard_thermal.h"
 %include "cahnhilliard_thermal_nodiffusion.h"
 %include "chparams.h"
 %include "run_ch_solver.hpp"
 %template(run_ch_vector_nonthermal) run_ch_solver<CHparamsVector , CahnHilliard2DRHS>; 
 %template(run_ch_scalar_nonthermal) run_ch_solver<CHparamsScalar , CahnHilliard2DRHS>;
 %template(run_ch_vector_thermal) run_ch_solver<CHparamsVector , CahnHilliard2DRHS_thermal>; 
 %template(run_ch_scalar_thermal) run_ch_solver<CHparamsScalar , CahnHilliard2DRHS_thermal>;
 %template(run_ch_vector_thermal_nodiffusion) run_ch_solver<CHparamsVector , CahnHilliard2DRHS_thermal_nodiffusion>; 
