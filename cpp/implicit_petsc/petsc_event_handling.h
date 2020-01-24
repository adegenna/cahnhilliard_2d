#ifndef __PETSC_EVENT_HANDLING_H__
#define __PETSC_EVENT_HANDLING_H__

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include "utils_ch_implicit.h"

void log_solution( Vec U , const std::string& outname );

PetscErrorCode EventFunction( TS ts , PetscReal t , Vec U , PetscScalar *fvalue , void *ctx );

PetscErrorCode PostEventFunction_ResetTemperatureGaussianProfile(TS ts,PetscInt nevents,PetscInt event_list[],PetscReal t,Vec U,PetscBool forwardsolve,void* ctx);

PetscErrorCode PostEventFunction_RecomputeThermalProperties(TS ts,PetscInt nevents,PetscInt event_list[],PetscReal t,Vec U,PetscBool forwardsolve,void* ctx);

void compute_new_temperature_profile( AppCtx* app , Vec U , PetscScalar T_amp , PetscScalar T_x , PetscScalar T_y , PetscScalar T_sigma  );

void compute_new_temperature_source_profile( AppCtx* app , Vec U , PetscScalar T_amp , PetscScalar T_x , PetscScalar T_y , PetscScalar T_sigma  );

#endif
