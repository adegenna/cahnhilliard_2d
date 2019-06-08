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

 %include "cahnhilliard.h"
