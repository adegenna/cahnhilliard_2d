#ifndef __CAHNHILLIARD_H__
#define __CAHNHILLIARD_H__

typedef std::vector<double> state_type;

struct ConstantValuedCahnHilliardParams
{
  double m;
  double gam;
  double b;
  double u;
  double alpha;
  double phi_star;
  double sigma;
};

struct VectorValuedCahnHilliardParams
{
  state_type m_xy;
  state_type gam_xy;
  state_type b_xy;
  state_type u_xy;
  state_type alpha_xy;
  state_type phi_star_xy;
  state_type sigma_xy;
};


class CahnHilliard2DRHS {

 public:

  CahnHilliard2DRHS(int nx, double dx);
  ~CahnHilliard2DRHS();
  void rhs_scalar_parameters(const state_type &c, state_type &dcdt, const double t);
  void rhs_field_parameters(const state_type &c, state_type &dcdt, const double t);
  void operator()(const state_type &c, state_type &dcdt, const double t);
  void setInitialConditions(state_type &x);
  double l2residual(const state_type&c);
  
 private:

  const int param_type_;
  
  double D_;
  double gamma_;
  double b_;
  double u_;
  double alpha_;
  double phi_star_;
  double sigma_;

  state_type D_xy_;
  state_type gamma_xy_;
  state_type b_xy_;
  state_type u_xy_;
  state_type alpha_xy_;
  state_type phi_star_xy_;
  state_type sigma_xy_;
  
  const int nx_;       // number of finite difference nodes in each dimension
  const double dx_;       // mesh size
  std::default_random_engine generator_;
  std::normal_distribution<double> noise_dist_;

  virtual const std::vector<double>& gammaCoefficient() const = 0;

  double laplace_component(double c);
  double laplace_component_field(const state_type& c, int i, int j);
  int idx2d_impl(int i, int j);
  int mod(int a, int b);
  int idx2d(int i, int j);
  
};


class ConstantCoefficientCahnHilliard : public CahnHilliard2DRHS {
 public:
  ConstantCoefficientCahnHilliard(const ConstantValuedCahnHilliardParams& p, int nx, double dx)
    : CahnHilliard2DRHS(nx, dx),
      gamma_(std::vector<double>(nx*nx, p.gam))
    {}

 private:
  const std::vector<double> gamma_;

  const std::vector<double>& gammaCoefficient() const override { return gamma_; }
};


class VariableCoefficientCahnHilliard : public CahnHilliard2DRHS {
 public:
  VariableCoefficientCahnHilliard(const VectorValuedCahnHilliardParams& p, int nx, double dx)
    : CahnHilliard2DRHS(nx, dx),
      params_(p)
    {}

 private:
  const VectorValuedCahnHilliardParams params_;

  const std::vector<double>& gammaCoefficient() const override { return params_.gam_xy; }
};




void write_state(const state_type &x , const int idx , const int nx );
void run_ch_solver_checkpointing(CHparams& chparams);
void run_ch_solver(CHparams& chparams);



#endif
