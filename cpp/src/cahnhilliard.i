 %module cahnhilliard

 %{
 /* Put header files here or function declarations like below */
 #include <iostream>
 #include <vector>
 #include <random>
 #include <fstream>
 #include <omp.h>
 #include "chparams.h"
 #include "run_ch_solver.h"
 %}

 %include <std_vector.i>
 %include typemaps.i
 %include std_string.i
 namespace std {
  %template(DoubleVector) vector<double>;
 };

 %include "chparams.h"
 %include "run_ch_solver.h"
