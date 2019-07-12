#ifndef __CHPARAMS_H__
#define __CHPARAMS_H__

#include <vector>
#include <algorithm>
#include <string>

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
  double dx, dy;
  int nx, ny;
  double dt_check;
  std::vector<double> x;
  std::string bc = "periodic";
  std::string rhs_type = "ch_non_thermal";
  double BC_dirichlet_ch;

  int idx2d(int i, int j);

 private:
  int idx2d_impl(int i, int j);
  int mod(int a, int b);
  
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
  double eps2_min, eps2_max, sigma_min, sigma_max, T_min, T_max, T_const;
  double sigma_noise;
  bool temperature_dependence = false;
  

  double compute_stability_limit(double dx , double dy);
  
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
  std::vector<double> T_const;
  double eps2_min, eps2_max, sigma_min, sigma_max, T_min, T_max;
  bool temperature_dependence = false;

  double sigma_noise;

  double compute_stability_limit(double dx , double dy);
  
};

#endif
