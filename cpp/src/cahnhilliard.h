#ifndef __CAHNHILLIARD_H__
#define __CAHNHILLIARD_H__

struct CHparams
{
  int param_type;
  
  double D;
  double gamma;
  double b;
  double u;
  double alpha;
  double phi_star;
  double sigma;

  std::vector<double> D_xy;
  std::vector<double> gamma_xy;
  std::vector<double> b_xy;
  std::vector<double> u_xy;
  std::vector<double> alpha_xy;
  std::vector<double> phi_star_xy;
  std::vector<double> sigma_xy;
  
  double t0;
  double tf;
  int iter = 0;
  double dx;
  double nx;
  double dt_check;
};

class CahnHilliard2DRHS {

 public:

  CahnHilliard2DRHS(CHparams& chp);
  ~CahnHilliard2DRHS();
  void rhs_scalar_parameters(const std::vector<double> &c, std::vector<double> &dcdt, const double t);
  void rhs_field_parameters(const std::vector<double> &c, std::vector<double> &dcdt, const double t);
  void operator()(const std::vector<double> &c, std::vector<double> &dcdt, const double t);
  void setInitialConditions(std::vector<double> &x);
  double l2residual(const std::vector<double> &c);
  
 private:

  CHparams& chp_;

  std::default_random_engine generator_;
  std::normal_distribution<double> noise_dist_;

  double laplace_component(double c);
  double laplace_component_field(const std::vector<double>& c, int i, int j);
  int idx2d_impl(int i, int j);
  int mod(int a, int b);
  int idx2d(int i, int j);
  
};

void write_state(const std::vector<double> &x , const int idx , const int nx );
void run_ch_solver_checkpointing(CHparams& chparams);
void run_ch_solver(CHparams& chparams);



#endif
