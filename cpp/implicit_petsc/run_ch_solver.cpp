#include <boost/numeric/odeint.hpp>
#include "stochastic_euler.hpp"
#include "cahnhilliard.h"
#include "cahnhilliard_thermal.h"
#include "cahnhilliard_thermal_nodiffusion.h"
#include "run_ch_solver.h"

void run_ch_solver_non_thermal( CHparamsVector& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS rhs = CahnHilliard2DRHS( chparams , info );
  
  std::vector<double> x;
  if (info.t0 == 0) {
    rhs.setInitialConditions(x);
    int iter = 0;
  }
  else {
    x        = info.x;
    int iter = info.iter;
  }
  
  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  const double stability_limit = chparams.compute_stability_limit(info.dx , info.dy); // just an estimate
  const double res0            = rhs.l2residual(x);

  std::cout << "residual at initial condition: " << res0 << std::endl;
  if (info.iter == 0)
    rhs.write_state(x,0,info.nx,info.ny,info.outdir);

  if (chparams.sigma_noise < 1e-2) {
    std::cout << "Solving deterministic (noise-free) CH" << std::endl;
    integrate_adaptive(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
    //boost::numeric::odeint::integrate_const(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
  }
  else {
    std::cout << "Solving stochastic CH" << std::endl;
    boost::mt19937 rng;
    boost::numeric::odeint::integrate_const( stochastic_euler() ,
		     std::make_pair( rhs , ornstein_stoch( rng , chparams.sigma_noise ) ),
		     x , info.t0 , info.tf , stability_limit/40. );
  }
  info.iter += 1;
  std::cout << "iter: " << info.iter << " , t = " << info.tf << ", relative residual: " << rhs.l2residual(x) / res0 << std::endl;
  rhs.write_state(x,info.iter,info.nx,info.ny,info.outdir);
  info.x = x;

};


void run_ch_solver_thermal_no_diffusion( CHparamsVector& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS_thermal_nodiffusion rhs = CahnHilliard2DRHS_thermal_nodiffusion( chparams , info );
  
  std::vector<double> x;
  if (info.t0 == 0) {
    rhs.setInitialConditions(x);
    int iter = 0;
  }
  else {
    x        = info.x;
    int iter = info.iter;
  }
  
  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  const double stability_limit = chparams.compute_stability_limit(info.dx , info.dy); // just an estimate
  const double res0            = rhs.l2residual(x);

  std::cout << "residual at initial condition: " << res0 << std::endl;
  if (info.iter == 0)
    rhs.write_state(x,0,info.nx,info.ny,info.outdir);

  if (chparams.sigma_noise < 1e-2) {
    std::cout << "Solving deterministic (noise-free) CH" << std::endl;
    integrate_adaptive(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
    //boost::numeric::odeint::integrate_const(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
  }
  else {
    std::cout << "Solving stochastic CH" << std::endl;
    boost::mt19937 rng;
    boost::numeric::odeint::integrate_const( stochastic_euler() ,
		     std::make_pair( rhs , ornstein_stoch( rng , chparams.sigma_noise ) ),
		     x , info.t0 , info.tf , stability_limit/40. );
  }
  info.iter += 1;
  std::cout << "iter: " << info.iter << " , t = " << info.tf << ", relative residual: " << rhs.l2residual(x) / res0 << std::endl;
  rhs.write_state(x,info.iter,info.nx,info.ny,info.outdir);
  info.x = x;

};


void run_ch_solver_thermal_with_diffusion( CHparamsVector& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS_thermal rhs = CahnHilliard2DRHS_thermal( chparams , info );
  
  std::vector<double> x;
  if (info.t0 == 0) {
    rhs.setInitialConditions(x);
    int iter = 0;
  }
  else {
    x        = info.x;
    int iter = info.iter;
  }
  
  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  const double stability_limit = chparams.compute_stability_limit(info.dx , info.dy); // just an estimate
  const double res0            = rhs.l2residual(x);

  std::cout << "residual at initial condition: " << res0 << std::endl;
  if (info.iter == 0)
    rhs.write_state(x,0,info.nx,info.ny,info.outdir);

  if (chparams.sigma_noise < 1e-2) {
    std::cout << "Solving deterministic (noise-free) CH" << std::endl;
    integrate_adaptive(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
    //boost::numeric::odeint::integrate_const(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
  }
  else {
    std::cout << "Solving stochastic CH" << std::endl;
    boost::mt19937 rng;
    boost::numeric::odeint::integrate_const( stochastic_euler() ,
		     std::make_pair( rhs , ornstein_stoch( rng , chparams.sigma_noise ) ),
		     x , info.t0 , info.tf , stability_limit/40. );
  }
  info.iter += 1;
  std::cout << "iter: " << info.iter << " , t = " << info.tf << ", relative residual: " << rhs.l2residual(x) / res0 << std::endl;
  rhs.write_state(x,info.iter,info.nx,info.ny,info.outdir);
  info.x = x;

};


void run_ch_solver_non_thermal( CHparamsScalar& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS rhs = CahnHilliard2DRHS( chparams , info );
  
  std::vector<double> x;
  if (info.t0 == 0) {
    rhs.setInitialConditions(x);
    int iter = 0;
  }
  else {
    x        = info.x;
    int iter = info.iter;
  }
  
  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  const double stability_limit = chparams.compute_stability_limit(info.dx , info.dy); // just an estimate
  const double res0            = rhs.l2residual(x);

  std::cout << "residual at initial condition: " << res0 << std::endl;
  if (info.iter == 0)
    rhs.write_state(x,0,info.nx,info.ny,info.outdir);

  if (chparams.sigma_noise < 1e-2) {
    std::cout << "Solving deterministic (noise-free) CH" << std::endl;
    integrate_adaptive(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
    //boost::numeric::odeint::integrate_const(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
  }
  else {
    std::cout << "Solving stochastic CH" << std::endl;
    boost::mt19937 rng;
    boost::numeric::odeint::integrate_const( stochastic_euler() ,
		     std::make_pair( rhs , ornstein_stoch( rng , chparams.sigma_noise ) ),
		     x , info.t0 , info.tf , stability_limit/40. );
  }
  info.iter += 1;
  std::cout << "iter: " << info.iter << " , t = " << info.tf << ", relative residual: " << rhs.l2residual(x) / res0 << std::endl;
  rhs.write_state(x,info.iter,info.nx,info.ny,info.outdir);
  info.x = x;

};


void run_ch_solver_thermal_no_diffusion( CHparamsScalar& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS_thermal_nodiffusion rhs = CahnHilliard2DRHS_thermal_nodiffusion( chparams , info );
  
  std::vector<double> x;
  if (info.t0 == 0) {
    rhs.setInitialConditions(x);
    int iter = 0;
  }
  else {
    x        = info.x;
    int iter = info.iter;
  }
  
  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  const double stability_limit = chparams.compute_stability_limit(info.dx , info.dy); // just an estimate
  const double res0            = rhs.l2residual(x);

  std::cout << "residual at initial condition: " << res0 << std::endl;
  if (info.iter == 0)
    rhs.write_state(x,0,info.nx,info.ny,info.outdir);

  if (chparams.sigma_noise < 1e-2) {
    std::cout << "Solving deterministic (noise-free) CH" << std::endl;
    integrate_adaptive(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
    //boost::numeric::odeint::integrate_const(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
  }
  else {
    std::cout << "Solving stochastic CH" << std::endl;
    boost::mt19937 rng;
    boost::numeric::odeint::integrate_const( stochastic_euler() ,
		     std::make_pair( rhs , ornstein_stoch( rng , chparams.sigma_noise ) ),
		     x , info.t0 , info.tf , stability_limit/40. );
  }
  info.iter += 1;
  std::cout << "iter: " << info.iter << " , t = " << info.tf << ", relative residual: " << rhs.l2residual(x) / res0 << std::endl;
  rhs.write_state(x,info.iter,info.nx,info.ny,info.outdir);
  info.x = x;

};


void run_ch_solver_thermal_with_diffusion( CHparamsScalar& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS_thermal rhs = CahnHilliard2DRHS_thermal( chparams , info );
  
  std::vector<double> x;
  if (info.t0 == 0) {
    rhs.setInitialConditions(x);
    int iter = 0;
  }
  else {
    x        = info.x;
    int iter = info.iter;
  }
  
  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  const double stability_limit = chparams.compute_stability_limit(info.dx , info.dy); // just an estimate
  const double res0            = rhs.l2residual(x);

  std::cout << "residual at initial condition: " << res0 << std::endl;
  if (info.iter == 0)
    rhs.write_state(x,0,info.nx,info.ny,info.outdir);

  if (chparams.sigma_noise < 1e-2) {
    std::cout << "Solving deterministic (noise-free) CH" << std::endl;
    integrate_adaptive(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
    //boost::numeric::odeint::integrate_const(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
  }
  else {
    std::cout << "Solving stochastic CH" << std::endl;
    boost::mt19937 rng;
    boost::numeric::odeint::integrate_const( stochastic_euler() ,
		     std::make_pair( rhs , ornstein_stoch( rng , chparams.sigma_noise ) ),
		     x , info.t0 , info.tf , stability_limit/40. );
  }
  info.iter += 1;
  std::cout << "iter: " << info.iter << " , t = " << info.tf << ", relative residual: " << rhs.l2residual(x) / res0 << std::endl;
  rhs.write_state(x,info.iter,info.nx,info.ny,info.outdir);
  info.x = x;

};

void run_ch_solver( CHparamsVector& chparams , SimInfo& info )
{
  
  if      (info.rhs_type.compare("ch_non_thermal") == 0) {
    run_ch_solver_non_thermal( chparams , info );
  }
  else if (info.rhs_type.compare("ch_thermal_no_diffusion") == 0) {
    run_ch_solver_thermal_no_diffusion( chparams , info );
  }
  else if (info.rhs_type.compare("ch_thermal_with_diffusion") == 0) {
    run_ch_solver_thermal_with_diffusion( chparams , info );
  }

};

void run_ch_solver( CHparamsScalar& chparams , SimInfo& info )
{
  
  if      (info.rhs_type.compare("ch_non_thermal") == 0) {
    run_ch_solver_non_thermal( chparams , info );
  }
  else if (info.rhs_type.compare("ch_thermal_no_diffusion") == 0) {
    run_ch_solver_thermal_no_diffusion( chparams , info );
  }
  else if (info.rhs_type.compare("ch_thermal_with_diffusion") == 0) {
    run_ch_solver_thermal_with_diffusion( chparams , info );
  }

};


void run_ch_solver_implicit( CHparamsVector& chparams , SimInfo& info ) {

  // Instantiate rhs
  CahnHilliard2DRHS rhs = CahnHilliard2DRHS( chparams , info );

  // Set initial conditions
  std::vector<double> x;
  if (info.t0 == 0) {
    rhs.setInitialConditions(x);
    int iter = 0;
  }
  else {
    x        = info.x;
    int iter = info.iter;
  }
  
  TS             ts;                   /* nonlinear solver */
  Vec            u,r;                  /* solution, residual vectors */
  Mat            J,Jmf = NULL;         /* Jacobian matrices */
  PetscErrorCode ierr;
  DM             da = NULL;
  SNES           snes;
  
  // *************** TEMPORARY HARD-CODED OPTIONS ***************
  double dt = 0.005;
  int Jtype = 2; /* Jacobian type
                    0: user provide Jacobian;
                    1: slow finite difference;
                    2: fd with coloring; */
  bool viewJacobian = PETSC_FALSE;
  double t_final    = 1.0;
  // ************************************************************
  
  PetscInitialize( NULL , NULL , "petsc_config.dat" , help );
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create distributed array (DMDA) to manage parallel grid and vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMDACreate2d(PETSC_COMM_WORLD, 
               DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,      // type of ghost nodes
               DMDA_STENCIL_BOX,                        // type of stencil
               11,11,                                   // global dimns of array
               PETSC_DECIDE,PETSC_DECIDE,               // #procs in each dimn
               1,                                       // DOF per node
               2,                                       // Stencil width
               NULL,NULL,&da);
  DMSetFromOptions(da);
  DMSetUp(da);
  
  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Extract global vectors from DMDA;
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMCreateGlobalVector(da,&u);
  VecDuplicate(u,&r);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create timestepping solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  TSCreate(PETSC_COMM_WORLD,&ts);
  TSSetProblemType(ts,TS_NONLINEAR);
  TSSetDM(ts,da);
  TSSetIFunction( ts , r , rhs );
  TSSetMaxTime( ts , t_final );
  TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set initial conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  FormInitialSolution(u,&user);
  TSSetSolution(ts,u);
  TSSetTimeStep(ts,dt);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set Jacobian evaluation routine
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr  = DMSetMatType(da,MATAIJ);
  ierr  = DMCreateMatrix(da,&J);
  TSGetSNES(ts,&snes);
  MatCreateSNESMF(snes,&Jmf);
  PetscPrintf(PETSC_COMM_WORLD,"JACOBIAN: matrix-free with coloring for finite-differences\n");
  SNESSetJacobian(snes,Jmf,J,SNESComputeJacobianDefaultColor,0);

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Sets various TS parameters from user options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  TSSetFromOptions(ts);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  TSSolve(ts,u);
  printf("Simulation done, press enter to continue...\n");
  getchar();
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  MatDestroy(&J);
  MatDestroy(&Jmf);
  VecDestroy(&u);
  VecDestroy(&r);
  TSDestroy(&ts);
  DMDestroy(&da);

  PetscFinalize();
  return ierr;

};
