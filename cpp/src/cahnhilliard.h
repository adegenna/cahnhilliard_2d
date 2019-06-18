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

  double D;
  double gamma;
  double b;
  double u;
  double alpha;
  double phi_star;
  double sigma;

  double compute_stability_limit(double dx)
  {
    return 0.5 * dx * dx * dx * dx / D / gamma;
  };
  
};

class CHparamsVector
{
  
 public:

  CHparamsVector()  { };
  ~CHparamsVector() { };
  
  std::vector<double> D;
  std::vector<double> gamma;
  std::vector<double> b;
  std::vector<double> u;
  std::vector<double> alpha;
  std::vector<double> phi_star;
  double sigma;

  double compute_stability_limit(double dx)
  {
    int idx_dmax = std::distance( D.begin()     , std::max_element( D.begin()     , D.end()     , abs_compare ) );
    int idx_gmax = std::distance( gamma.begin() , std::max_element( gamma.begin() , gamma.end() , abs_compare ) );
    double dmax = D[     idx_dmax ];
    double gmax = gamma[ idx_gmax ];
    return 0.5 * dx * dx * dx * dx / dmax / gmax;
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
			   const std::vector<double>& D ,
			   const std::vector<double>& u ,
			   const std::vector<double>& b);
  int idx2d_impl(int i, int j);
  int mod(int a, int b);
  int idx2d(int i, int j);
  
};

void write_state(const std::vector<double> &x , const int idx , const int nx );




#endif
