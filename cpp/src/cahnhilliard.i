 %module cahnhilliard
 %{
 /* Put header files here or function declarations like below */
 #include <iostream>
 #include <vector>
 #include <random>
 #include <fstream>
 #include <omp.h>
 #include <boost/numeric/odeint.hpp>
 #include "cahnhilliard.h"
 %}

 %include <std_vector.i>
 %include typemaps.i
 namespace std {
  %template(DoubleVector) vector<double>;
 };
 %include "cahnhilliard.h"
