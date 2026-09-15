program validate_ml_surrogate
  !
  ! Standalone bit-verification of the ML surrogate forward pass.  Loads the
  ! trained weights directly via the HDF5 Fortran API (no PFLOTRAN init
  ! plumbing), then evaluates every row of the golden test-vector CSV and
  ! compares against the expected columns.  Both output stages are handled,
  ! selected by the weight file (root attribute output_stage_id):
  !   stage 0 (sinh rate):  expected_rate_per_yr
  !   stage 1 (node_min_v1): expected_d_vf_csh, expected_d_s_csh, expected_d_vf_m4ac,
  !            expected_d_vf_other, expected_dphi, expected_cs_new, expected_sink_{Si,Ca,Mg,Al,HCO3,TotalH}
  ! Tolerance: relative 1e-9 with an absolute floor of 1e-10 (a cell exactly at
  ! the residual porosity has a reference rate of 0 while round-off can give
  ! ~1e-9; see NEURAL_ODE_HANDOFF.md section 6).
  !
  ! Run:
  !   make validate_ml_surrogate
  !   ./validate_ml_surrogate <weights.h5> <test_vectors.csv>
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24 (both output stages, header-driven columns: 2026-09-14)
  !

#include "petsc/finclude/petscsys.h"
  use petscsys
  use hdf5
  use ML_Surrogate_Aux_module

  implicit none

  integer, parameter :: NEXP = 12
  character(len=24), parameter :: exp_names(NEXP) = [character(len=24) :: &
    'expected_d_vf_csh', 'expected_d_s_csh', 'expected_d_vf_m4ac', 'expected_d_vf_other', &
    'expected_dphi', 'expected_cs_new', 'expected_sink_Si', 'expected_sink_Ca', &
    'expected_sink_Mg', 'expected_sink_Al', 'expected_sink_HCO3', 'expected_sink_TotalH']

  type(ml_mlp_type) :: mlp
  type(ml_node_min_out_type) :: nm
  character(len=1024) :: filename, csvname
  character(len=8192) :: line, hdr
  character(len=512) :: errmsg
  PetscReal :: rate, expected, vals(64), rel, worst, got(NEXP)
  PetscReal, allocatable :: x_in(:)
  PetscErrorCode :: ierr
  integer :: nargs, hdf5_err, ios, i, k, nrows, ncol, erate_col, exp_col(NEXP), lerr
  integer(HID_T) :: file_id

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRA(ierr)

  nargs = command_argument_count()
  if (nargs < 2) then
    write(*,'(a)') 'Usage: validate_ml_surrogate <weights.h5> <test_vectors.csv>'
    call PetscFinalize(ierr); CHKERRA(ierr)
    stop 1
  endif
  call get_command_argument(1, filename)
  call get_command_argument(2, csvname)

  call h5open_f(hdf5_err)
  call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, hdf5_err)
  if (hdf5_err < 0) then
    write(*,'(a,a,a)') 'Cannot open HDF5 file "', trim(filename), '".'
    call PetscFinalize(ierr); CHKERRA(ierr)
    stop 2
  endif

  call MLPCreate(mlp)
  call MLPLoadFromFileId(mlp, file_id, lerr, errmsg)
  if (lerr /= 0) then
    write(*,'(a,a)') 'Weight file rejected: ', trim(errmsg)
    call PetscFinalize(ierr); CHKERRA(ierr)
    stop 5
  endif
  write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') 'loaded: n_in=', mlp%n_in, ' hidden=', mlp%n_h0, '/', mlp%n_h1, '/', &
    mlp%n_h2, ' n_out=', mlp%n_out, ' output_stage=', mlp%stage

  call h5fclose_f(file_id, hdf5_err)
  call h5close_f(hdf5_err)

  open(unit=10, file=trim(csvname), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    write(*,'(a,a,a)') 'Cannot open CSV "', trim(csvname), '".'
    call PetscFinalize(ierr); CHKERRA(ierr)
    stop 3
  endif

  ! Header-driven column layout: features are cols 4..3+n_in (after eval_id,
  ! cell_x_m, time_y); the expected columns are located by name.
  read(10,'(a)',iostat=ios) hdr
  ncol = 1
  do i = 1, len_trim(hdr)
    if (hdr(i:i) == ',') ncol = ncol + 1
  enddo
  erate_col = ColumnOf(hdr, 'expected_rate_per_yr')
  do k = 1, NEXP
    exp_col(k) = ColumnOf(hdr, trim(exp_names(k)))
  enddo
  if (mlp%stage == ML_STAGE_PHI_SINH .and. erate_col == 0) then
    write(*,'(a)') 'CSV has no expected_rate_per_yr column (stage 0 file).'
    stop 6
  endif
  if (mlp%stage == ML_STAGE_NODE_MIN .and. any(exp_col == 0)) then
    write(*,'(a)') 'CSV lacks one of the expected_* columns of node_min_v1.'
    stop 6
  endif
  allocate(x_in(mlp%n_in))

  worst = 0.d0
  nrows = 0
  do
    read(10,'(a)',iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) == 0) cycle
    do i = 1, len_trim(line)
      if (line(i:i) == ',') line(i:i) = ' '
    enddo
    read(line,*,iostat=ios) vals(1:ncol)
    if (ios /= 0) cycle

    x_in = vals(4:3+mlp%n_in)        ! features, idx order 1..n_in
    select case(mlp%stage)
      case(ML_STAGE_PHI_SINH)
        expected = vals(erate_col)
        call MLPForward(mlp, x_in, rate)
        rel = abs(rate - expected) / max(abs(expected), 1.d-30)
        worst = max(worst, rel)
      case(ML_STAGE_NODE_MIN)
        call MLPForwardNodeMin(mlp, x_in, nm)
        got = [nm%d_vf_csh, nm%d_s_csh, nm%d_vf_m4ac, nm%d_vf_other, nm%dphi, nm%cs_new, &
               nm%sink_Si, nm%sink_Ca, nm%sink_Mg, nm%sink_Al, nm%sink_HCO3, nm%sink_TotalH]
        do k = 1, NEXP
          expected = vals(exp_col(k))
          rel = abs(got(k) - expected) / max(abs(expected), 1.d-10)
          worst = max(worst, rel)
        enddo
    end select
    nrows = nrows + 1
  enddo
  close(10)

  write(*,'(i0,a,es12.4)') nrows, ' vectors checked; worst relative error = ', worst
  if (worst < 1.d-9) then
    write(*,'(a)') 'PASS -- Fortran forward pass matches the golden vectors.'
  else
    write(*,'(a)') 'FAIL -- transform chain broken (check log10 mask, clip '// &
                   'bounds, feature order, TotalH sign, output-stage attributes).'
  endif

  call MLPDestroy(mlp)
  call PetscFinalize(ierr); CHKERRA(ierr)
  if (worst >= 1.d-9) stop 4

contains

! ************************************************************************** !

function ColumnOf(hdr, name) result(col)
  !
  ! 1-based index of the comma-separated header field equal to name; 0 if absent.
  !
  ! Author: Claude Fable 5.1 (for Heeho Park)
  ! Date: 2026-09-14
  !
  implicit none
  character(len=*) :: hdr, name
  integer :: col
  integer :: i, start, n

  col = 0; n = 1; start = 1
  do i = 1, len_trim(hdr) + 1
    if (i > len_trim(hdr) .or. hdr(i:i) == ',') then
      if (hdr(start:i-1) == name) then
        col = n; return
      endif
      n = n + 1; start = i + 1
    endif
  enddo
end function ColumnOf

end program validate_ml_surrogate
