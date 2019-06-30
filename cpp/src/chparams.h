#ifndef __CHPARAMS_H__
#define __CHPARAMS_H__

#include <vector>
#include <algorithm>

static bool abs_compare(int a, int b)
{
  return (std::abs(a) < std::abs(b));
}

class SimInfo
{
  
 public:

  SimInfo() { };
  ~SimInfo() { };
  
  double t0;
  double tf;
  int iter = 0;
  double dx;
  int nx;
  double dt_check;
  std::vector<double> x;

};

class CHparamsScalar
{
  
 public:

  CHparamsScalar()  { };
  ~CHparamsScalar() { };

  double eps_2;
  double b;
  double u;
  double sigma;
  double m;
  double DT;
  double f_T;
  double eps2_min, eps2_max, sigma_min, sigma_max, T_min, T_max;
  double sigma_noise;
  bool temperature_dependence = false;

  double compute_stability_limit(double dx);
  
};

class CHparamsVector
{
  
 public:

  CHparamsVector()  { };
  ~CHparamsVector() { };
  
  std::vector<double> eps_2;
  std::vector<double> b;
  std::vector<double> u;
  std::vector<double> sigma;
  std::vector<double> m;
  std::vector<double> DT;
  std::vector<double> f_T;
  double eps2_min, eps2_max, sigma_min, sigma_max, T_min, T_max;
  bool temperature_dependence = false;

  double sigma_noise;

  double compute_stability_limit(double dx);
  
};

#endif
