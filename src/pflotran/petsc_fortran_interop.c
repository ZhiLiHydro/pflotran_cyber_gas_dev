/*
  Fortran bindings for the TS driver (pmc_ts_driver.F90, which declares
  their interfaces) of TSSetEventHandler and SNESSetComputeInitialGuess:
  PETSc has no C stubs for them (3.24.5 through main, which declares only
  their interfaces). Names are mangled as PETSc's own stubs are.

  C trampolines call the Fortran callbacks with every argument by
  reference (a PETSc object's Fortran handle holds only its C pointer).
  The TS/SNES owns each wrapper context through a composed container.

  The post-event callback returns in dt_post the first step after the
  event (> 0), passed to TSSetPostEventStep so that TSAdapt may reject
  it; with PETSc's default (PETSC_DECIDE) it never does.

  Author: Hansol Suh
  Date: 09/22/26
*/
#include <petscts.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
  #define pflotran_ts_set_event_handler_           PFLOTRAN_TS_SET_EVENT_HANDLER
  #define pflotran_snes_set_compute_initial_guess_ PFLOTRAN_SNES_SET_COMPUTE_INITIAL_GUESS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
  #define pflotran_ts_set_event_handler_           pflotran_ts_set_event_handler
  #define pflotran_snes_set_compute_initial_guess_ pflotran_snes_set_compute_initial_guess
#endif

typedef void (*PFLOTRANEventIndicatorFn)(TS *, PetscReal *, Vec *, PetscReal *, void *, PetscErrorCode *);
typedef void (*PFLOTRANPostEventFn)(TS *, PetscInt *, PetscInt *, PetscReal *, Vec *, PetscBool *, PetscReal *, void *, PetscErrorCode *);

typedef struct {
  PFLOTRANEventIndicatorFn indicator;
  PFLOTRANPostEventFn      postevent;
  void                    *fctx;
} PFLOTRANEventCtx;

static PetscErrorCode PFLOTRANEventIndicator(TS ts, PetscReal t, Vec U, PetscReal fvalue[], void *ctx)
{
  PFLOTRANEventCtx *c    = (PFLOTRANEventCtx *)ctx;
  PetscErrorCode    ierr = PETSC_SUCCESS;

  PetscFunctionBegin;
  c->indicator(&ts, &t, &U, fvalue, c->fctx, &ierr);
  PetscCall(ierr);
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PFLOTRANPostEvent(TS ts, PetscInt nevents_zero, PetscInt events_zero[], PetscReal t, Vec U, PetscBool forwardsolve, void *ctx)
{
  PFLOTRANEventCtx *c       = (PFLOTRANEventCtx *)ctx;
  PetscErrorCode    ierr    = PETSC_SUCCESS;
  PetscReal         dt_post = -1.0;
  PetscInt          dummy   = 0;

  PetscFunctionBegin;
  c->postevent(&ts, &nevents_zero, nevents_zero > 0 ? events_zero : &dummy, &t, &U, &forwardsolve, &dt_post, c->fctx, &ierr);
  PetscCall(ierr);
  if (dt_post > 0.0) PetscCall(TSSetPostEventStep(ts, dt_post));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PFLOTRANEventCtxDestroy(void **ctx)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(*ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PETSC_EXTERN void pflotran_ts_set_event_handler_(TS *ts, PetscInt *nevents, PetscInt direction[], PetscBool terminate[], PFLOTRANEventIndicatorFn indicator, PFLOTRANPostEventFn postevent, void *fctx, PetscErrorCode *ierr)
{
  PFLOTRANEventCtx *c;
  PetscContainer    container;

  *ierr = PetscNew(&c);
  if (*ierr) return;
  c->indicator = indicator;
  c->postevent = postevent;
  c->fctx      = fctx;
  *ierr = PetscContainerCreate(PetscObjectComm((PetscObject)*ts), &container);
  if (*ierr) return;
  *ierr = PetscContainerSetPointer(container, c);
  if (*ierr) return;
  *ierr = PetscContainerSetCtxDestroy(container, PFLOTRANEventCtxDestroy);
  if (*ierr) return;
  *ierr = PetscObjectCompose((PetscObject)*ts, "PFLOTRANEventCtx", (PetscObject)container);
  if (*ierr) return;
  *ierr = PetscContainerDestroy(&container);
  if (*ierr) return;
  *ierr = TSSetEventHandler(*ts, *nevents, direction, terminate, PFLOTRANEventIndicator, PFLOTRANPostEvent, c);
}

/*
  SNESSetComputeInitialGuess: the Fortran callback receives the SNES, the
  initial guess (modifiable) and its context.
*/
typedef void (*PFLOTRANInitialGuessFn)(SNES *, Vec *, void *, PetscErrorCode *);

typedef struct {
  PFLOTRANInitialGuessFn func;
  void                  *fctx;
} PFLOTRANInitialGuessCtx;

static PetscErrorCode PFLOTRANInitialGuess(SNES snes, Vec X, void *ctx)
{
  PFLOTRANInitialGuessCtx *c    = (PFLOTRANInitialGuessCtx *)ctx;
  PetscErrorCode           ierr = PETSC_SUCCESS;

  PetscFunctionBegin;
  c->func(&snes, &X, c->fctx, &ierr);
  PetscCall(ierr);
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PFLOTRANInitialGuessCtxDestroy(void **ctx)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(*ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PETSC_EXTERN void pflotran_snes_set_compute_initial_guess_(SNES *snes, PFLOTRANInitialGuessFn func, void *fctx, PetscErrorCode *ierr)
{
  PFLOTRANInitialGuessCtx *c;
  PetscContainer           container;

  *ierr = PetscNew(&c);
  if (*ierr) return;
  c->func = func;
  c->fctx = fctx;
  *ierr   = PetscContainerCreate(PetscObjectComm((PetscObject)*snes), &container);
  if (*ierr) return;
  *ierr = PetscContainerSetPointer(container, c);
  if (*ierr) return;
  *ierr = PetscContainerSetCtxDestroy(container, PFLOTRANInitialGuessCtxDestroy);
  if (*ierr) return;
  *ierr = PetscObjectCompose((PetscObject)*snes, "PFLOTRANInitialGuessCtx", (PetscObject)container);
  if (*ierr) return;
  *ierr = PetscContainerDestroy(&container);
  if (*ierr) return;
  *ierr = SNESSetComputeInitialGuess(*snes, PFLOTRANInitialGuess, c);
}
