#ifndef __PETSC_EVENT_HANDLING_H__
#define __PETSC_EVENT_HANDLING_H__

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>

PetscErrorCode EventFunction( TS ts , PetscReal t , Vec U , PetscScalar *fvalue , void *ctx );

PetscErrorCode PostEventFunction(TS ts,PetscInt nevents,PetscInt event_list[],PetscReal t,Vec U,PetscBool forwardsolve,void* ctx);

#endif
