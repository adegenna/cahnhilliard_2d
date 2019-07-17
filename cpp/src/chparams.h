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
  
  double t0, tf;                // Initial and final simulation times
  int iter = 0;                 // Current simulation iteration
  double dx, dy;                // Spatial discretization
  int nx, ny;                   // Number of grid points in (x,y)
  std::vector<double> x;        // Optional, used to specify current state IFF t0 != 0 (MAKE THIS MORE ROBUST)
  std::string bc = "periodic";  // Boundary condition type: "periodic", "neumann", "dirichlet", "mixed_neumann_bottom_dirichlet", "mixed_neumann_top_dirichlet"
  std::string rhs_type = "ch_non_thermal"; // RHS type: "ch_non_thermal", "ch_thermal_no_diffusion", "ch_thermal_with_diffusion"
  double BC_dirichlet_ch;       // Used to specify BC value for dirichlet BC

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
  double L_kuhn, N, L_omega, X_min, X_max; // Polymer parameters
  double sigma_noise = 0.0;
  bool temperature_dependence = false;

  double compute_stability_limit(double dx , double dy);
  double convert_temperature_to_flory_huggins( const double T ,
					       const double T_min ,
					       const double T_max ,
					       const double X_min ,
					       const double X_max );
  double compute_eps2_from_polymer_params(     const double T ,
					       const double m ,
					       const double L_kuhn ,
					       const double N );
  double compute_sigma_from_polymer_params(    const double T ,
					       const double m ,
					       const double L_kuhn ,
					       const double L_omega ,
					       const double N );
  void compute_and_set_eps2_and_sigma_from_polymer_params( const double T );

};

class CHparamsVector
{
  
 public:

  CHparamsVector() { };
  CHparamsVector(int nx , int ny);
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
  double L_kuhn, N, L_omega, X_min, X_max; // Polymer parameters
  bool temperature_dependence = false;
  double sigma_noise = 0.0;

  double compute_stability_limit(double dx , double dy);
  double convert_temperature_to_flory_huggins( const double T ,
					       const double T_min ,
					       const double T_max ,
					       const double X_min ,
					       const double X_max );
  double compute_eps2_from_polymer_params(     const double T ,
					       const double m ,
					       const double L_kuhn ,
					       const double N );
  double compute_sigma_from_polymer_params(    const double T ,
					       const double m ,
					       const double L_kuhn ,
					       const double L_omega ,
					       const double N );
  void compute_and_set_eps2_and_sigma_from_polymer_params( const double T ,
							   SimInfo& info );
  
  
};

#endif
