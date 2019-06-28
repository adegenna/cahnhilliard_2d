#ifndef __CAHNHILLIARD_H__
#define __CAHNHILLIARD_H__

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
  double sigma_noise;

  double compute_stability_limit(double dx)
  {
    return 0.5 * dx * dx * dx * dx / eps_2;
  };
  
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
  double sigma_noise;

  double compute_stability_limit(double dx)
  {
    int idx_gmax = std::distance( eps_2.begin() , std::max_element( eps_2.begin() , eps_2.end() , abs_compare ) );
    double gmax = eps_2[ idx_gmax ];
    return 0.5 * dx * dx * dx * dx / gmax;
  };

  
};

class CahnHilliard2DRHS {

 public:

  CahnHilliard2DRHS(CHparamsScalar& chp , SimInfo& info);
  CahnHilliard2DRHS(CHparamsVector& chp , SimInfo& info);
  ~CahnHilliard2DRHS();
  void rhs(const std::vector<double> &c, std::vector<double> &dcdt, const double t);
  void operator()(const std::vector<double> &c, std::vector<double> &dcdt, const double t);
  void setInitialConditions(std::vector<double> &x);
  double l2residual(const std::vector<double> &c);
  
 private:

  CHparamsVector chpV_;
  SimInfo& info_;

  std::default_random_engine generator_;
  std::normal_distribution<double> noise_dist_;
  
  double laplace_component(int i ,
			   const std::vector<double>& c ,
			   const std::vector<double>& u ,
			   const std::vector<double>& b);
  int idx2d_impl(int i, int j);
  int mod(int a, int b);
  int idx2d(int i, int j);
  
};

void write_state(const std::vector<double> &x , const int idx , const int nx );




#endif
