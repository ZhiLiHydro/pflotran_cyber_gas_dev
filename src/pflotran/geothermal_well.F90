module Geothermal_Well_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Well_class
  use PM_Base_class
  use Option_module

  implicit none ! TODO: look this up

  private 

  ! Geothermal Quasi-Implicit Well Model
  type, public, extends(pm_well_qi_type) :: pm_well_geoth_qi_type
    contains
    
    !   procedure, public :: ReadPMBlock => PMWellReadPMBlockGeoThQI
    !   procedure, public :: ReadSimulationOptionsBlock => &
    !           PMWellReadSimOptionsBlockGeoThQI
    !   procedure, public :: Setup => PMWellSetupGeoThQI
    !   procedure, public :: Solve => PMWellSolveGeoThQI
    !   procedure, public :: SolveFlow => PMWellSolveFlowGeoThQI
    !   procedure, public :: ModifyFlowResidual => PMWellModifyFlowResGeoThQI
    !   procedure, public :: ModifyFlowJacobian => PMWellModifyFlowJacGeoThQI
    !   procedure, public :: UpdateFlowRates => PMWellUpdateFlowRatesGeoThQI
    !   procedure, public :: InitializeTimestep => PMWellInitializeTimestepGeoThQI
    !   procedure, public :: UpdateFlowProperties => UpdateFlowPropertiesGeoThQI
  end type pm_well_geoth_qi_type

!   public :: PMWellGeoThQICreate
  
  contains

end module Geothermal_Well_class