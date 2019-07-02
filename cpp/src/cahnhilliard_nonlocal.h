#ifndef __CAHNHILLIARD_NONLOCAL_H__
#define __CAHNHILLIARD_NONLOCAL_H__

#include <vector>
#include "chparams.h"

void compute_ch_nonlocal(const std::vector<double> &c,
			 std::vector<double> &dcdt,
			 const double t,
			 CHparamsVector chpV,
			 SimInfo& info);

double laplace_component(int i ,
                         const std::vector<double>& c ,
                         const std::vector<double>& u ,
                         const std::vector<double>& b );


#endif
