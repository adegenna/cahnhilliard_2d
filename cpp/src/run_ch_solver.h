#include "chparams.h"

void run_ch_solver( CHparamsVector& chparams , SimInfo& info );
void run_ch_solver( CHparamsScalar& chparams , SimInfo& info );

void run_ch_solver_non_thermal( CHparamsVector& chparams , SimInfo& info );
void run_ch_solver_thermal_no_diffusion( CHparamsVector& chparams , SimInfo& info );
void run_ch_solver_thermal_with_diffusion( CHparamsVector& chparams , SimInfo& info );

void run_ch_solver_non_thermal( CHparamsScalar& chparams , SimInfo& info );
void run_ch_solver_thermal_no_diffusion( CHparamsScalar& chparams , SimInfo& info );
void run_ch_solver_thermal_with_diffusion( CHparamsScalar& chparams , SimInfo& info );
