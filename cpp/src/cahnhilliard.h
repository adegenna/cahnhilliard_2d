#ifndef __CAHNHILLIARD_H__
#define __CAHNHILLIARD_H__

typedef std::vector<double> state_type;

struct CHparams
{
  double m;
  double gam;
  double b;
  double u;
  double alpha;
  double phi_star;
  double sigma;
};

class CahnHilliard2DRHS {

 public:

  CahnHilliard2DRHS(CHparams& chp, int nx, double dx);
  ~CahnHilliard2DRHS();
  void operator()(const state_type &c, state_type &dcdt, const double t);
  void setInitialConditions(state_type &x);
  double l2residual(const state_type&c);
  
 private:

  const double D_;     // diffusion coefficient
  const double gamma_; // the other term
  const double b_;
  const double u_;
  const double alpha_;
  const double phi_star_;
  const int nx_;       // number of finite difference nodes in each dimension
  const double dx_;       // mesh size
  const double sigma_;
  std::default_random_engine generator_;
  std::normal_distribution<double> noise_dist_;

  double laplace_component(double c);
  int idx2d_impl(int i, int j);
  int mod(int a, int b);
  int idx2d(int i, int j);
  
};

void write_state(const state_type &x , const int idx , const int nx );
void run_ch_solver(CHparams& chparams, const int nx, const double dx, const int checkpoint, const int maxsteps);



#endif
