module TS_Optimizer_class
  !
  ! Calibration of parameters with PETSc TAO: minimizes an objective whose
  ! value and gradient an evaluator returns (for PMC_TS_Driver_class, a
  ! rerun of the forward run and its adjoint sweep). Physics-free: it sees
  ! parameter values only.
  !
  ! The optimization variables are the parameters, or their log10
  ! (TRANSFORM LOG10, the default for most that must be positive);
  ! bounds are given in the parameters' units. A point where a parameter
  ! that must be positive is not (or 10**q is not a positive double) is
  ! not evaluated. An evaluation without a gradient returns +Inf: Armijo's
  ! line search, the default here in place of More-Thuente, backtracks
  ! from it, More-Thuente stops. An evaluation that asks to stop
  ! (TS_OPT_EVAL_STOP, e.g. at the wall-clock limit) ends the optimization
  ! at TAO's last accepted point: later evaluations return +Inf without
  ! running. The last
  ! evaluation is kept, so TAO's first point, the forward run's own, is
  ! not evaluated again. Every evaluation is a row of
  ! <prefix>-optimization.log, with the parameters (the screen shows the
  ! objective only). TAO's options take the prefix adj_
  ! (-adj_tao_monitor, -adj_tao_test_gradient, ...).
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
#include "petsc/finclude/petsctao.h"
  use petsctao

  use PFLOTRAN_Constants_module
  use Option_module

  implicit none

  private

  PetscInt, parameter, public :: TS_OPT_TRANSFORM_DEFAULT = 0
  PetscInt, parameter, public :: TS_OPT_TRANSFORM_NONE = 1
  PetscInt, parameter, public :: TS_OPT_TRANSFORM_LOG10 = 2

  ! evaluator status: a value, or stop the optimization; others mean no
  ! value at that point
  PetscInt, parameter, public :: TS_OPT_EVAL_OK = 0
  PetscInt, parameter, public :: TS_OPT_EVAL_STOP = -1
  ! not evaluated: outside the parameters' domain
  PetscInt, parameter :: TS_OPT_EVAL_OUTSIDE = -2

  ! what returns the objective and its gradient at given parameters
  type, public, abstract :: ts_evaluator_type
  contains
    procedure(TSOptEvaluateIface), deferred :: Evaluate
  end type ts_evaluator_type

  abstract interface
    subroutine TSOptEvaluateIface(this,p,f,g,status)
      ! f, the objective at parameters p, and g = df/dp, on every process;
      ! status: TS_OPT_EVAL_OK, TS_OPT_EVAL_STOP, or another (no value).
      use petscsys
      import :: ts_evaluator_type
      implicit none
      class(ts_evaluator_type) :: this
      PetscReal :: p(:)
      PetscReal :: f
      PetscReal :: g(:)
      PetscInt :: status
    end subroutine TSOptEvaluateIface
  end interface

  ! what TAO hands the callbacks
  type :: ts_optimizer_ctx_type
    class(ts_optimizer_type), pointer :: optimizer => null()
  end type ts_optimizer_ctx_type

  type, public :: ts_optimizer_type
    type(option_type), pointer :: option => null()
    ! the OPTIMIZATION block: TAO type, iteration limit, gradient
    ! tolerances (TAO's defaults if uninitialized), and the settings of
    ! parameters by number
    character(len=MAXSTRINGLENGTH) :: tao_type = 'bqnls'
    PetscInt :: max_iterations = UNINITIALIZED_INTEGER
    PetscReal :: gatol = UNINITIALIZED_DOUBLE
    PetscReal :: grtol = UNINITIALIZED_DOUBLE
    PetscInt :: nspec = 0
    PetscInt, allocatable :: spec_param(:)
    PetscInt, allocatable :: spec_transform(:)
    PetscReal, allocatable :: spec_lower(:)
    PetscReal, allocatable :: spec_upper(:)
    ! the parameters (Setup): names, whether each must be positive,
    ! transforms, bounds and initial values
    PetscInt :: n = 0
    character(len=MAXSTRINGLENGTH), allocatable :: name(:)
    PetscBool, allocatable :: positive(:)
    PetscInt, allocatable :: transform(:)
    PetscReal, allocatable :: lower(:)
    PetscReal, allocatable :: upper(:)
    PetscReal, allocatable :: p0(:)
    ! the last evaluation: variables q, objective, gradient df/dp, status
    PetscBool :: have_last = PETSC_FALSE
    PetscReal, allocatable :: q_last(:)
    PetscReal :: f_last = 0.d0
    PetscReal, allocatable :: g_last(:)
    PetscInt :: status_last = 0
    PetscInt :: num_evaluations = 0
    PetscInt :: num_failures = 0
    ! an evaluation asked to stop
    PetscBool :: stopped = PETSC_FALSE
    ! the evaluation log (IO rank; newunit units are negative)
    PetscBool :: log_open = PETSC_FALSE
    PetscInt :: log_unit = 0
    ! the result
    PetscInt :: iterations = 0
    PetscInt :: reason = 0
    PetscReal :: f_opt = 0.d0
    PetscReal, allocatable :: p_opt(:)
    class(ts_evaluator_type), pointer :: evaluator => null()
    type(ts_optimizer_ctx_type), pointer :: ctx => null()
  contains
    procedure, public :: AddParameterSettings => TSOptAddParameterSettings
    procedure, public :: Setup => TSOptSetup
    procedure, public :: SetInitialEvaluation => TSOptSetInitialEvaluation
    procedure, public :: Solve => TSOptSolve
    procedure, public :: RegressionOutput => TSOptRegressionOutput
    procedure, public :: Destroy => TSOptDestroy
  end type ts_optimizer_type

  public :: TSOptimizerCreate

contains

! ************************************************************************** !

function TSOptimizerCreate()
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  type(ts_optimizer_type), pointer :: TSOptimizerCreate

  allocate(TSOptimizerCreate)
  allocate(TSOptimizerCreate%spec_param(0))
  allocate(TSOptimizerCreate%spec_transform(0))
  allocate(TSOptimizerCreate%spec_lower(0))
  allocate(TSOptimizerCreate%spec_upper(0))

end function TSOptimizerCreate

! ************************************************************************** !

subroutine TSOptAddParameterSettings(this,iparam,transform,lower,upper)
  !
  ! The transform and bounds of parameter iparam (a number of the
  ! parameter list), from the input; resolved at Setup.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_optimizer_type) :: this
  PetscInt :: iparam
  PetscInt :: transform
  PetscReal :: lower
  PetscReal :: upper

  this%nspec = this%nspec + 1
  this%spec_param = [this%spec_param, iparam]
  this%spec_transform = [this%spec_transform, transform]
  this%spec_lower = [this%spec_lower, lower]
  this%spec_upper = [this%spec_upper, upper]

end subroutine TSOptAddParameterSettings

! ************************************************************************** !

subroutine TSOptSetup(this,option,name,positive,log_scale,p0)
  !
  ! The parameters: names, whether each must be positive, whether LOG10
  ! is its default transform, and their initial values; resolves and
  ! checks the settings.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_optimizer_type) :: this
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: name(:)
  PetscBool :: positive(:)
  PetscBool :: log_scale(:)
  PetscReal :: p0(:)

  PetscInt :: i, j

  this%option => option
  this%n = size(p0)
  allocate(this%name(this%n))
  allocate(this%positive(this%n))
  allocate(this%transform(this%n))
  allocate(this%lower(this%n))
  allocate(this%upper(this%n))
  allocate(this%p0(this%n))
  allocate(this%q_last(this%n))
  allocate(this%g_last(this%n))
  allocate(this%p_opt(this%n))
  this%name = name
  this%positive = positive
  this%p0 = p0
  this%p_opt = p0
  this%lower = -huge(1.d0)
  this%upper = huge(1.d0)
  do j = 1, this%n
    if (log_scale(j)) then
      this%transform(j) = TS_OPT_TRANSFORM_LOG10
    else
      this%transform(j) = TS_OPT_TRANSFORM_NONE
    endif
  enddo
  do i = 1, this%nspec
    j = this%spec_param(i)
    if (j < 1 .or. j > this%n) then
      write(option%io_buffer,'(a,i0,a,i0,a)') 'ADJOINT,OPTIMIZATION: &
        &PARAMETER ',j,' does not exist (the ADJOINT block has ',this%n, &
        ' parameters).'
      call PrintErrMsg(option)
    endif
    if (any(this%spec_param(1:i-1) == j)) then
      write(option%io_buffer,'(a,i0,a)') 'ADJOINT,OPTIMIZATION: &
        &PARAMETER ',j,' is given twice.'
      call PrintErrMsg(option)
    endif
    if (this%spec_transform(i) /= TS_OPT_TRANSFORM_DEFAULT) &
      this%transform(j) = this%spec_transform(i)
    if (Initialized(this%spec_lower(i))) this%lower(j) = this%spec_lower(i)
    if (Initialized(this%spec_upper(i))) this%upper(j) = this%spec_upper(i)
  enddo
  do j = 1, this%n
    if (positive(j) .and. .not.(p0(j) > 0.d0)) then
      option%io_buffer = 'ADJOINT,OPTIMIZATION: parameter "' // &
        trim(name(j)) // '" must be positive.'
      call PrintErrMsg(option)
    endif
    if (this%transform(j) == TS_OPT_TRANSFORM_LOG10) then
      if (.not.(p0(j) > 0.d0)) then
        option%io_buffer = 'ADJOINT,OPTIMIZATION: parameter "' // &
          trim(name(j)) // '" is not positive; use TRANSFORM NONE.'
        call PrintErrMsg(option)
      endif
      this%lower(j) = max(this%lower(j),0.d0)
    endif
    if (.not.(this%lower(j) < this%upper(j)) .or. &
        p0(j) < this%lower(j) .or. p0(j) > this%upper(j)) then
      option%io_buffer = 'ADJOINT,OPTIMIZATION: the initial value of &
        &parameter "' // trim(name(j)) // '" is not inside its BOUNDS.'
      call PrintErrMsg(option)
    endif
  enddo

end subroutine TSOptSetup

! ************************************************************************** !

subroutine TSOptSetInitialEvaluation(this,f,g)
  !
  ! The forward run's own objective and gradient, at the initial
  ! parameters: TAO's first point.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

  implicit none

  class(ts_optimizer_type) :: this
  PetscReal :: f
  PetscReal :: g(:)

  if (.not.(ieee_is_finite(f) .and. all(ieee_is_finite(g)))) then
    call PrintErrMsg(this%option,'ADJOINT,OPTIMIZATION: the objective or &
      &its gradient is not finite at the initial parameters.')
  endif
  call TSOptVariables(this,this%p0,this%q_last)
  this%f_last = f
  this%g_last = g
  this%status_last = 0
  this%have_last = PETSC_TRUE
  this%num_evaluations = 1
  call TSOptLogOpen(this)
  call TSOptPrintEvaluation(this,this%p0)

end subroutine TSOptSetInitialEvaluation

! ************************************************************************** !

subroutine TSOptSolve(this,evaluator)
  !
  ! Minimizes the objective from the initial parameters; p_opt, f_opt.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_optimizer_type), target :: this
  class(ts_evaluator_type), target :: evaluator

  type(option_type), pointer :: option
  Tao :: tao
  TaoType :: tao_type
  TaoLineSearch :: linesearch
  character(len=MAXWORDLENGTH) :: ls_type
  PetscBool :: flag
  Vec :: x, xl, xu
  TaoConvergedReason :: reason
  PetscReal :: q(this%n), ql(this%n), qu(this%n), gnorm, cnorm, xdiff
  PetscReal :: gatol, grtol
  PetscInt :: nlocal, j
  PetscErrorCode :: ierr

  option => this%option
  this%evaluator => evaluator
  if (.not.associated(this%ctx)) allocate(this%ctx)
  this%ctx%optimizer => this

  ! the variables: all on the first process
  nlocal = 0
  if (option%myrank == 0) nlocal = this%n
  call VecCreateMPI(option%mycomm,nlocal,this%n,x,ierr);CHKERRQ(ierr)
  call VecDuplicate(x,xl,ierr);CHKERRQ(ierr)
  call VecDuplicate(x,xu,ierr);CHKERRQ(ierr)
  call TSOptVariables(this,this%p0,q)
  do j = 1, this%n
    ql(j) = PETSC_NINFINITY
    qu(j) = PETSC_INFINITY
    select case(this%transform(j))
      case(TS_OPT_TRANSFORM_LOG10)
        if (this%lower(j) > 0.d0) ql(j) = log10(this%lower(j))
        if (this%upper(j) < huge(1.d0)) qu(j) = log10(this%upper(j))
      case default
        if (this%lower(j) > -huge(1.d0)) ql(j) = this%lower(j)
        if (this%upper(j) < huge(1.d0)) qu(j) = this%upper(j)
    end select
  enddo
  call TSOptSetVec(this,x,q)
  call TSOptSetVec(this,xl,ql)
  call TSOptSetVec(this,xu,qu)

  call TaoCreate(option%mycomm,tao,ierr);CHKERRA(ierr)
  call TaoSetOptionsPrefix(tao,'adj_',ierr);CHKERRA(ierr)
  call TaoSetType(tao,trim(this%tao_type),ierr);CHKERRA(ierr)
  call TaoSetSolution(tao,x,ierr);CHKERRA(ierr)
  call TaoSetVariableBounds(tao,xl,xu,ierr);CHKERRA(ierr)
  call TaoSetObjectiveAndGradient(tao,PETSC_NULL_VEC, &
                                  TSOptObjectiveAndGradient,this%ctx, &
                                  ierr);CHKERRA(ierr)
  ! CHKERRA: CHKERRQ would return silently with a half-configured TAO
  if (Initialized(this%max_iterations)) then
    call TaoSetMaximumIterations(tao,this%max_iterations,ierr);CHKERRA(ierr)
  endif
  if (Initialized(this%gatol)) then
    gatol = this%gatol
    grtol = this%grtol
    call TaoSetTolerances(tao,gatol,grtol,PETSC_CURRENT_REAL, &
                          ierr);CHKERRA(ierr)
  endif
  call TaoSetFromOptions(tao,ierr);CHKERRA(ierr)
  ! -adj_tao_type may have changed it
  call TaoGetType(tao,tao_type,ierr);CHKERRA(ierr)
  this%tao_type = tao_type
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,'adj_','-tao_ls_type',flag, &
                           ierr);CHKERRA(ierr)
  if (.not.flag) then
    ! a method without a line search (e.g. nm) has none; SetType resets the
    ! line search's options
    call TaoGetLineSearch(tao,linesearch,ierr);CHKERRA(ierr)
    if (.not.PetscObjectIsNull(linesearch)) then
      call TaoLineSearchGetType(linesearch,ls_type,ierr);CHKERRA(ierr)
      if (trim(ls_type) == TAOLINESEARCHMT) then
        call TaoLineSearchSetType(linesearch,TAOLINESEARCHARMIJO, &
                                  ierr);CHKERRA(ierr)
        call TaoLineSearchSetFromOptions(linesearch,ierr);CHKERRA(ierr)
      endif
    endif
  endif

  write(option%io_buffer,'(a,i0,a)') ' OPTIMIZATION: TAO ' // &
    trim(this%tao_type) // ' over ',this%n,' parameter(s); evaluations in ' &
    // trim(option%global_prefix) // trim(option%group_prefix) // &
    '-optimization.log'
  call PrintMsg(option)
  call TaoSolve(tao,ierr);CHKERRA(ierr)

  call TaoGetSolutionStatus(tao,this%iterations,this%f_opt,gnorm,cnorm, &
                            xdiff,reason,ierr);CHKERRQ(ierr)
  this%reason = reason%v
  call TSOptGetVec(this,x,q)
  call TSOptParameters(this,q,this%p_opt)
  write(option%io_buffer,'(a,i0,a,i0,a,i0,a,a,a,es21.13,a,es12.5)') &
    ' OPTIMIZATION: ',this%iterations,' iterations, ', &
    this%num_evaluations,' evaluations (',this%num_failures,' failed); ', &
    trim(TSOptReasonName(this%reason)),'; f = ',this%f_opt,', |g| = ',gnorm
  call PrintMsg(option)
  do j = 1, this%n
    write(option%io_buffer,'(a,a,a,es21.13)') ' OPTIMIZATION: ', &
      trim(this%name(j)),' = ',this%p_opt(j)
    call PrintMsg(option)
  enddo
  if (this%stopped) then
    write(option%io_buffer,'(a,i0,a)') 'OPTIMIZATION stopped at the &
      &wall-clock limit (evaluation ',this%num_evaluations,'); the &
      &parameters above are TAO''s last accepted point.'
    call PrintWrnMsg(option)
  else if (this%reason < 0) then
    call PrintWrnMsg(option,'OPTIMIZATION: TAO did not converge; the &
      &output is that of its last point.')
  endif
  if (this%log_open) then
    close(this%log_unit)
    this%log_open = PETSC_FALSE
  endif

  call TaoDestroy(tao,ierr);CHKERRQ(ierr)
  call VecDestroy(x,ierr);CHKERRQ(ierr)
  call VecDestroy(xl,ierr);CHKERRQ(ierr)
  call VecDestroy(xu,ierr);CHKERRQ(ierr)
  nullify(this%evaluator)

end subroutine TSOptSolve

! ************************************************************************** !

subroutine TSOptObjectiveAndGradient(tao,x,f,g,ctx,ierr)
  !
  ! TAO's objective and gradient: f and df/dq at the variables x.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  Tao :: tao
  Vec :: x
  PetscReal :: f
  Vec :: g
  type(ts_optimizer_ctx_type) :: ctx
  PetscErrorCode :: ierr

  class(ts_optimizer_type), pointer :: this
  PetscReal, allocatable :: q(:), gq(:)

  this => ctx%optimizer
  allocate(q(this%n),gq(this%n))
  call TSOptGetVec(this,x,q)
  call TSOptEvaluate(this,q,f,gq)
  call TSOptSetVec(this,g,gq)
  ierr = 0

end subroutine TSOptObjectiveAndGradient

! ************************************************************************** !

subroutine TSOptEvaluate(this,q,f,gq)
  !
  ! The objective and df/dq at variables q (every process), from the last
  ! evaluation if q is its point; +Inf and zero without a value.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_optimizer_type) :: this
  PetscReal :: q(:)
  PetscReal :: f
  PetscReal :: gq(:)

  PetscReal :: p(this%n), dp_dq(this%n)

  call TSOptParameters(this,q,p,dp_dq)
  if (this%stopped) then
    f = PETSC_INFINITY
    gq = 0.d0
    return
  endif
  if (.not.(this%have_last .and. all(q == this%q_last))) then
    if (TSOptInDomain(this,q)) then
      call this%evaluator%Evaluate(p,this%f_last,this%g_last, &
                                   this%status_last)
    else
      this%status_last = TS_OPT_EVAL_OUTSIDE
    endif
    this%q_last = q
    this%have_last = PETSC_TRUE
    this%num_evaluations = this%num_evaluations + 1
    if (this%status_last == TS_OPT_EVAL_STOP) then
      this%stopped = PETSC_TRUE
    else if (this%status_last /= TS_OPT_EVAL_OK) then
      this%num_failures = this%num_failures + 1
    endif
    call TSOptPrintEvaluation(this,p)
  endif
  if (this%status_last /= TS_OPT_EVAL_OK) then
    f = PETSC_INFINITY
    gq = 0.d0
  else
    f = this%f_last
    gq = this%g_last*dp_dq
  endif

end subroutine TSOptEvaluate

! ************************************************************************** !

subroutine TSOptPrintEvaluation(this,p)
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_optimizer_type) :: this
  PetscReal :: p(:)

  call TSOptLogEvaluation(this,p)
  if (this%status_last == TS_OPT_EVAL_STOP) then
    write(this%option%io_buffer,'(a,i5,a)') &
      ' OPTIMIZATION evaluation ',this%num_evaluations, &
      ': stopped at the wall-clock limit'
  else if (this%status_last == TS_OPT_EVAL_OUTSIDE) then
    write(this%option%io_buffer,'(a,i5,a)') &
      ' OPTIMIZATION evaluation ',this%num_evaluations, &
      ': not run, a parameter is outside its domain'
  else if (this%status_last /= TS_OPT_EVAL_OK) then
    write(this%option%io_buffer,'(a,i5,a,i0,a)') &
      ' OPTIMIZATION evaluation ',this%num_evaluations, &
      ': no gradient (status ',this%status_last,')'
  else
    write(this%option%io_buffer,'(a,i5,a,es21.13,a,es12.5)') &
      ' OPTIMIZATION evaluation ',this%num_evaluations,': f = ', &
      this%f_last,', |df/dp| = ',sqrt(sum(this%g_last**2))
  endif
  call PrintMsg(this%option)

end subroutine TSOptPrintEvaluation

! ************************************************************************** !

subroutine TSOptLogOpen(this)
  !
  ! <prefix>-optimization.log: one row per evaluation, its number, status
  ! (0: a value, -1: stopped, -2: not run, outside the parameters'
  ! domain, else none), objective, gradient norm, parameters and gradient
  ! (NaN without a value).
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_optimizer_type) :: this

  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: j

  if (.not.OptionIsIORank(this%option)) return
  filename = trim(this%option%global_prefix) // &
             trim(this%option%group_prefix) // '-optimization.log'
  open(newunit=this%log_unit,file=filename,action='write',status='replace')
  this%log_open = PETSC_TRUE
  write(this%log_unit,'(a)') '# ADJOINT OPTIMIZATION evaluations: number, &
    &status (0: value, -1: stopped, -2: not run, outside the domain, else &
    &none), f, |df/dp|, p_j, df/dp_j'
  do j = 1, this%n
    write(this%log_unit,'(a,i0,a,a)') '# p_',j,': ',trim(this%name(j))
  enddo

end subroutine TSOptLogOpen

! ************************************************************************** !

subroutine TSOptLogEvaluation(this,p)
  !
  ! The last evaluation's row of the log.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

  implicit none

  class(ts_optimizer_type) :: this
  PetscReal :: p(:)

  PetscReal :: f, gnorm, g(this%n)

  if (.not.this%log_open) return
  if (this%status_last == TS_OPT_EVAL_OK) then
    f = this%f_last
    g = this%g_last
    gnorm = sqrt(sum(g**2))
  else
    f = ieee_value(f,ieee_quiet_nan)
    g = f
    gnorm = f
  endif
  write(this%log_unit,'(i6,1x,i4,*(1x,es24.16))') this%num_evaluations, &
    this%status_last,f,gnorm,p,g
  flush(this%log_unit)

end subroutine TSOptLogEvaluation

! ************************************************************************** !

subroutine TSOptVariables(this,p,q)
  !
  ! The variables q of parameters p.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_optimizer_type) :: this
  PetscReal :: p(:)
  PetscReal :: q(:)

  PetscInt :: j

  do j = 1, this%n
    if (this%transform(j) == TS_OPT_TRANSFORM_LOG10) then
      q(j) = log10(p(j))
    else
      q(j) = p(j)
    endif
  enddo

end subroutine TSOptVariables

! ************************************************************************** !

subroutine TSOptParameters(this,q,p,dp_dq)
  !
  ! The parameters p of variables q, and dp/dq; the initial parameters
  ! exactly at their own variables.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_optimizer_type) :: this
  PetscReal :: q(:)
  PetscReal :: p(:)
  PetscReal, optional :: dp_dq(:)

  PetscReal :: q0(this%n)
  PetscInt :: j

  call TSOptVariables(this,this%p0,q0)
  do j = 1, this%n
    if (this%transform(j) == TS_OPT_TRANSFORM_LOG10) then
      ! 0 and +Inf outside the doubles (TSOptInDomain)
      if (q(j) > log10(huge(1.d0))) then
        p(j) = PETSC_INFINITY
      else if (q(j) < log10(tiny(1.d0))) then
        p(j) = 0.d0
      else
        p(j) = 10.d0**q(j)
      endif
      if (q(j) == q0(j)) p(j) = this%p0(j)
      if (present(dp_dq)) dp_dq(j) = log(10.d0)*p(j)
    else
      p(j) = q(j)
      if (present(dp_dq)) dp_dq(j) = 1.d0
    endif
  enddo

end subroutine TSOptParameters

! ************************************************************************** !

function TSOptInDomain(this,q)
  !
  ! Whether variables q are parameters that can be evaluated: finite, a
  ! parameter that must be positive is, and 10**q is a positive double.
  !
  ! Author: Hansol Suh
  ! Date: 10/11/26
  !
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

  implicit none

  class(ts_optimizer_type) :: this
  PetscReal :: q(:)
  PetscBool :: TSOptInDomain

  PetscInt :: j

  TSOptInDomain = all(ieee_is_finite(q))
  do j = 1, this%n
    if (.not.TSOptInDomain) exit
    if (this%transform(j) == TS_OPT_TRANSFORM_LOG10) then
      TSOptInDomain = q(j) >= log10(tiny(1.d0)) .and. &
                      q(j) <= log10(huge(1.d0))
    else if (this%positive(j)) then
      TSOptInDomain = q(j) > 0.d0
    endif
  enddo

end function TSOptInDomain

! ************************************************************************** !

subroutine TSOptGetVec(this,v,q)
  !
  ! q = v on every process (v's entries are on the first).
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_optimizer_type) :: this
  Vec :: v
  PetscReal :: q(:)

  PetscReal, pointer :: v_p(:)
  PetscErrorCode :: ierr

  call VecGetArrayRead(v,v_p,ierr);CHKERRQ(ierr)
  if (this%option%myrank == 0) q = v_p(1:this%n)
  call VecRestoreArrayRead(v,v_p,ierr);CHKERRQ(ierr)
  call MPI_Bcast(q,this%n,MPI_DOUBLE_PRECISION,ZERO_INTEGER_MPI, &
                 this%option%mycomm,ierr);CHKERRQ(ierr)

end subroutine TSOptGetVec

! ************************************************************************** !

subroutine TSOptSetVec(this,v,q)
  !
  ! v = q, the same on every process.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_optimizer_type) :: this
  Vec :: v
  PetscReal :: q(:)

  PetscReal, pointer :: v_p(:)
  PetscErrorCode :: ierr

  call VecGetArray(v,v_p,ierr);CHKERRQ(ierr)
  if (this%option%myrank == 0) v_p(1:this%n) = q
  call VecRestoreArray(v,v_p,ierr);CHKERRQ(ierr)

end subroutine TSOptSetVec

! ************************************************************************** !

function TSOptReasonName(reason)
  !
  ! TaoConvergedReason as text.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  PetscInt :: reason
  character(len=MAXSTRINGLENGTH) :: TSOptReasonName

  select case(reason)
    case(3)
      TSOptReasonName = 'converged: gradient absolute tolerance'
    case(4)
      TSOptReasonName = 'converged: gradient relative tolerance'
    case(5)
      TSOptReasonName = 'converged: gradient reduction'
    case(6)
      TSOptReasonName = 'converged: step tolerance'
    case(7)
      TSOptReasonName = 'converged: objective minimum'
    case(-2)
      TSOptReasonName = 'diverged: maximum iterations'
    case(-4)
      TSOptReasonName = 'diverged: Inf or NaN'
    case(-5)
      TSOptReasonName = 'diverged: maximum function evaluations'
    case(-6)
      TSOptReasonName = 'diverged: line search failure'
    case(-7)
      TSOptReasonName = 'diverged: trust region reduction'
    case default
      write(TSOptReasonName,'(a,i0)') 'TaoConvergedReason ',reason
  end select

end function TSOptReasonName

! ************************************************************************** !

subroutine TSOptRegressionOutput(this,fid_out)
  !
  ! Writes the result as a section of the regression file.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_optimizer_type) :: this
  PetscInt :: fid_out

  PetscInt :: j

  if (.not.OptionIsIORank(this%option)) return
  write(fid_out,'(''-- GENERIC: Optimization --'')')
  write(fid_out,'(3x,a,'': '',es21.13)') 'iterations', &
    dble(this%iterations)
  write(fid_out,'(3x,a,'': '',es21.13)') 'converged reason', dble(this%reason)
  write(fid_out,'(3x,a,'': '',es21.13)') 'f', this%f_opt
  do j = 1, this%n
    write(fid_out,'(3x,a,i0,a,'': '',es21.13)') 'p_',j,'*', this%p_opt(j)
  enddo

end subroutine TSOptRegressionOutput

! ************************************************************************** !

subroutine TSOptDestroy(this)
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_optimizer_type) :: this

  if (this%log_open) then
    close(this%log_unit)
    this%log_open = PETSC_FALSE
  endif
  if (associated(this%ctx)) then
    nullify(this%ctx%optimizer)
    deallocate(this%ctx)
    nullify(this%ctx)
  endif
  nullify(this%evaluator)

end subroutine TSOptDestroy

end module TS_Optimizer_class
