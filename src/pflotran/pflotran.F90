!=======================================================================
! Please see LICENSE and COPYRIGHT files at top of repository.
!=======================================================================

! sjb - sleep fn for MPI debugging
module Fortran_Sleep
  use, intrinsic :: iso_c_binding, only: c_int
  implicit none

  interface
     !  should be unsigned int ... not available in Fortran
     !  OK until highest bit gets set.
     function FortSleep (seconds)  bind ( C, name="sleep" )
         import
         integer (c_int) :: FortSleep
         integer (c_int), intent (in), VALUE :: seconds
     end function FortSleep
  end interface
end module Fortran_Sleep




program pflotran

  use Driver_class
  use Simulation_Base_class
  use Factory_PFLOTRAN_module
  use PFLOTRAN_Constants_module

  use Fortran_Sleep

  implicit none

#include "petsc/finclude/petscsys.h"

  class(simulation_base_type), pointer :: simulation
  class(driver_type), pointer :: driver
  PetscInt :: iflag

  ! For MPI debugging, set i!=0 here and then run with MPI.
  ! Then attach debugger instances to each running process
  ! and modify value of i to 0 with debugger to proceed to dbg
  integer :: i, duration
  i = 0
  do while (i.ne.0)
    write(*, "('  Sleeping... ',i4)") i
    duration = FortSleep(i)
  end do


  driver => DriverCreate()
  call FactoryPFLOTRANInitialize(driver,simulation)
  call simulation%InitializeRun()
  if (driver%status == PROCEED) then
    call simulation%ExecuteRun()
  endif
  call simulation%FinalizeRun()
  call SimulationBaseDestroy(simulation)
  call FactoryPFLOTRANFinalize(driver)
  iflag = driver%exit_code
  call DriverDestroy(driver)
  call exit(iflag)

end program pflotran
