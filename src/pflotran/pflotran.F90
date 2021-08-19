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
  
  use Option_module
  use Simulation_Base_class
  use Multi_Simulation_module
  use Factory_PFLOTRAN_module
  use Factory_Subsurface_module
  use PFLOTRAN_Constants_module
  use PFLOTRAN_Provenance_module, only : PrintProvenanceToScreen
  
  use Fortran_Sleep

  implicit none

#include "petsc/finclude/petscsys.h"

  class(simulation_base_type), pointer :: simulation
  ! multisimulation enables multiple simulations to be run concurrently
  ! and/or one after another until a specified set of simulations has 
  ! completed.
  type(multi_simulation_type), pointer :: multisimulation
  type(option_type), pointer :: option
  
  integer :: i, duration
  
  ! For MPI debugging, set i!=0 here and then run with MPI.
  ! Then attach debugger instances to each running process
  ! and modify value of i to 0 with debugger to proceed to dbg
  i = 0
  do while (i.ne.0)
    write(*, "('  Sleeping... ',i4)") i
    duration = FortSleep(i)
  end do


  nullify(simulation)
  nullify(multisimulation)
  option => OptionCreate()
  call OptionInitMPI(option)
  call PFLOTRANInitializePrePetsc(multisimulation,option)
  call OptionInitPetsc(option)
  if (option%myrank == option%io_rank .and. option%print_to_screen) then
    !call PrintProvenanceToScreen()
  endif

  do ! multi-simulation loop
    call PFLOTRANInitializePostPetsc(simulation,multisimulation,option)

    call simulation%InitializeRun()

    if (option%status == PROCEED) then
      call simulation%ExecuteRun()
    endif

    call simulation%FinalizeRun()
    call simulation%Strip()
    deallocate(simulation)
    nullify(simulation)
    call PFLOTRANFinalize(option)
    if (MultiSimulationDone(multisimulation)) exit
  enddo ! multi-simulation loop
  call OptionFinalize(option)

end program pflotran
