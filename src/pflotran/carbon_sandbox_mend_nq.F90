module Carbon_Sandbox_MEND_nQ_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Carbon_Sandbox_Base_class
  use Carbon_Sandbox_MEND_class
  use Option_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  ! The 2013 MEND model of Wang et al. (doi:10.1890/12-0681.1) carries a single
  ! sorbed carbon pool Q that exchanges with DOC through adsorption (F6) and
  ! desorption (F7). This sandbox partitions Q into n mineral-specific pools
  ! Q_1 ... Q_n so that soil mineralogy controls sorption capacity and rate,
  ! after Iyer, MacFarlane, and Zavarin. Only F6 and F7 change; every other
  ! flux and pool equation is identical to Carbon_Sandbox_MEND_class.
  !
  ! Per pool i:
  !   F6_i = k_ads,i*(1 - Q_i/Q_max,i)*D
  !   F7_i = k_des,i*(Q_i/Q_max,i)
  !   dQ_i = F6_i - F7_i
  ! and the DOC balance sums over pools: dD = ... - (sum(F6_i) - sum(F7_i))
  !
  ! Per-pool parameters are either given directly (Q_MAX, K_ADS, K_DES) or
  ! derived from mineral properties. For a derived pool:
  !   Q_max,i = w_i*SA_i*SD_i*MW_OM*f_C*1.d21/N_A        [mg C/g soil]
  !   k_ads,i = K_ads_base*Q_max,i/sum_j(Q_max,j)
  !   K_eq,i  = K_d,i*(1-phi)*w_i /
  !             (phi*S_w*sum_j(w_j/rho_j)*Q_max,i)       [g soil/mg C]
  !   k_des,i = k_ads,i/K_eq,i
  ! K_eq,i depends on porosity and water saturation, which PFLOTRAN tracks as
  ! live state, so k_des,i is evaluated per cell rather than once at setup.

  ! atoms per mole. Defined locally, as in pm_ert.F90 and pm_ufd_biosphere.F90.
  PetscReal, parameter :: AVOGADRO_NUMBER = 6.02214076d23

  type, public :: carbon_sandbox_q_pool_type
    character(len=MAXWORDLENGTH) :: site_name
    character(len=MAXWORDLENGTH) :: species_name
    ! supplied directly
    PetscReal :: Q_max                 ! mg C/g soil
    PetscReal :: K_ads                 ! kg/kg-s
    PetscReal :: K_des                 ! kg/kg-s
    ! supplied as mineral properties
    PetscReal :: mass_fraction         ! unitless
    PetscReal :: specific_surface_area ! m^2/g
    PetscReal :: site_density          ! sites/nm^2
    PetscReal :: Kd_lscie              ! mL/g
    PetscReal :: mineral_density       ! g/cm^3
    type(carbon_sandbox_q_pool_type), pointer :: next
  end type carbon_sandbox_q_pool_type

  type, extends(carbon_sandbox_mend_type), public :: carbon_sandbox_mend_nq_type
    ! populated while reading; flattened into the arrays below at setup
    type(carbon_sandbox_q_pool_type), pointer :: q_pool_list
    PetscInt :: n_q_pool
    PetscInt, pointer :: Q_pool_index(:)
    PetscReal, pointer :: Q_pool_max(:)
    PetscReal, pointer :: Q_pool_kads(:)
    PetscReal, pointer :: Q_pool_kdes(:)
    PetscReal, pointer :: Q_pool_Keq_coef(:)
    PetscBool, pointer :: Q_pool_derived(:)
    ! mineralogy parameters shared by all derived pools
    PetscReal :: MW_OM                 ! g/mol
    PetscReal :: f_C                   ! unitless
    PetscReal :: K_ads_base            ! kg/kg-s
  contains
    procedure, public :: ReadInput => CarbonMENDnQReadInput
    procedure, public :: Setup => CarbonMENDnQSetup
    procedure, public :: Evaluate => CarbonMENDnQEvaluate
    procedure, public :: Strip => CarbonMENDnQStrip
  end type carbon_sandbox_mend_nq_type

  public :: CarbonMENDnQCreate

contains

! ************************************************************************** !

function CarbonMENDnQCreate()
  !
  ! Allocates and initializes the carbon sandbox object
  !
  ! Author: Richard Tran Mills
  ! Date: 10/26/26
  !
  class(carbon_sandbox_mend_nq_type), pointer :: CarbonMENDnQCreate

  class(carbon_sandbox_mend_nq_type), pointer :: this

  allocate(this)
  call CarbonMENDInit(this)
  ! the single-Q members inherited from the MEND sandbox are unused here
  nullify(this%q_pool_list)
  this%n_q_pool = 0
  nullify(this%Q_pool_index)
  nullify(this%Q_pool_max)
  nullify(this%Q_pool_kads)
  nullify(this%Q_pool_kdes)
  nullify(this%Q_pool_Keq_coef)
  nullify(this%Q_pool_derived)
  this%MW_OM = UNINITIALIZED_DOUBLE
  this%f_C = UNINITIALIZED_DOUBLE
  this%K_ads_base = UNINITIALIZED_DOUBLE

  CarbonMENDnQCreate => this

end function CarbonMENDnQCreate

! ************************************************************************** !

function CarbonQPoolCreate()
  !
  ! Allocates and initializes a single sorbed carbon pool
  !
  ! Author: Richard Tran Mills
  ! Date: 10/26/26
  !
  type(carbon_sandbox_q_pool_type), pointer :: CarbonQPoolCreate

  type(carbon_sandbox_q_pool_type), pointer :: this

  allocate(this)
  this%site_name = ''
  this%species_name = ''
  this%Q_max = UNINITIALIZED_DOUBLE
  this%K_ads = UNINITIALIZED_DOUBLE
  this%K_des = UNINITIALIZED_DOUBLE
  this%mass_fraction = UNINITIALIZED_DOUBLE
  this%specific_surface_area = UNINITIALIZED_DOUBLE
  this%site_density = UNINITIALIZED_DOUBLE
  this%Kd_lscie = UNINITIALIZED_DOUBLE
  this%mineral_density = UNINITIALIZED_DOUBLE
  nullify(this%next)

  CarbonQPoolCreate => this

end function CarbonQPoolCreate

! ************************************************************************** !

subroutine CarbonMENDnQReadInput(this,input,option)
  !
  ! Reads parameters from input deck
  !
  ! Author: Richard Tran Mills
  ! Date: 10/26/26
  !
  use Input_Aux_module
  use Option_module
  use String_module

  class(carbon_sandbox_mend_nq_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: err_string
  PetscBool :: found

  err_string = 'CHEMISTRY,CARBON_SANDBOX,MEND_NQ'
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',err_string)
    call StringToUpper(keyword)

    found = PETSC_FALSE
    call CarbonBaseReadSelectCase(this,input,keyword,found,err_string,option)
    if (found) cycle

    call CarbonMENDReadSelectCase(this,input,keyword,found,err_string,option)
    if (found) cycle

    select case(keyword)
      case('MINERALOGY')
        call CarbonMENDnQReadMineralogy(this,input,option)
      case('SORPTION_SITE')
        call CarbonMENDnQReadSorptionSite(this,input,option)
      case default
        call InputKeywordUnrecognized(input,keyword,err_string,option)
    end select
  enddo

end subroutine CarbonMENDnQReadInput

! ************************************************************************** !

subroutine CarbonMENDnQReadMineralogy(this,input,option)
  !
  ! Reads the mineralogy parameters shared by all derived sorption sites
  !
  ! Author: Richard Tran Mills
  ! Date: 10/26/26
  !
  use Input_Aux_module
  use Option_module
  use String_module

  class(carbon_sandbox_mend_nq_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXWORDLENGTH) :: internal_units
  character(len=MAXSTRINGLENGTH) :: err_string

  err_string = 'CHEMISTRY,CARBON_SANDBOX,MEND_NQ,MINERALOGY'

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',err_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
      case('MW_OM')
        call InputReadDouble(input,option,this%MW_OM)
        call InputErrorMsg(input,option,keyword,err_string)
        ! the Q_max expression of Iyer et al. is posed in g/mol
        internal_units = 'g/mol'
        call InputReadAndConvertUnits(input,this%MW_OM,internal_units, &
                                      keyword,option)
      case('F_C')
        call InputReadDouble(input,option,this%f_C)
        call InputErrorMsg(input,option,keyword,err_string)
      case('K_ADS_BASE')
        call InputReadDouble(input,option,this%K_ads_base)
        call InputErrorMsg(input,option,keyword,err_string)
        internal_units = 'kg/kg-s'
        call InputReadAndConvertUnits(input,this%K_ads_base,internal_units, &
                                      keyword,option)
      case default
        call InputKeywordUnrecognized(input,keyword,err_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine CarbonMENDnQReadMineralogy

! ************************************************************************** !

subroutine CarbonMENDnQReadSorptionSite(this,input,option)
  !
  ! Reads a single sorbed carbon pool and appends it to the list
  !
  ! Author: Richard Tran Mills
  ! Date: 10/26/26
  !
  use Input_Aux_module
  use Option_module
  use String_module

  class(carbon_sandbox_mend_nq_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXWORDLENGTH) :: internal_units
  character(len=MAXSTRINGLENGTH) :: err_string
  type(carbon_sandbox_q_pool_type), pointer :: new_pool, cur_pool

  new_pool => CarbonQPoolCreate()
  call InputReadWord(input,option,new_pool%site_name,PETSC_TRUE)
  call InputErrorMsg(input,option,'name', &
                     'CHEMISTRY,CARBON_SANDBOX,MEND_NQ,SORPTION_SITE')
  err_string = 'CHEMISTRY,CARBON_SANDBOX,MEND_NQ,SORPTION_SITE,' // &
               trim(new_pool%site_name)

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',err_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
      case('Q_SPECIES_NAME')
        call InputReadWord(input,option,new_pool%species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('Q_MAX')
        call InputReadDouble(input,option,new_pool%Q_max)
        call InputErrorMsg(input,option,keyword,err_string)
        ! compared directly against the state variables, which
        ! CarbonMENDnQEvaluate carries in mg C/g soil
        internal_units = 'mg/g'
        call InputReadAndConvertUnits(input,new_pool%Q_max,internal_units, &
                                      keyword,option)
      case('K_ADS')
        call InputReadDouble(input,option,new_pool%K_ads)
        call InputErrorMsg(input,option,keyword,err_string)
        internal_units = 'kg/kg-s'
        call InputReadAndConvertUnits(input,new_pool%K_ads,internal_units, &
                                      keyword,option)
      case('K_DES')
        call InputReadDouble(input,option,new_pool%K_des)
        call InputErrorMsg(input,option,keyword,err_string)
        internal_units = 'kg/kg-s'
        call InputReadAndConvertUnits(input,new_pool%K_des,internal_units, &
                                      keyword,option)
      case('MASS_FRACTION')
        call InputReadDouble(input,option,new_pool%mass_fraction)
        call InputErrorMsg(input,option,keyword,err_string)
      case('SITE_DENSITY') ! sites/nm^2
        call InputReadDouble(input,option,new_pool%site_density)
        call InputErrorMsg(input,option,keyword,err_string)
      case('SPECIFIC_SURFACE_AREA')
        call InputReadDouble(input,option, &
                             new_pool%specific_surface_area)
        call InputErrorMsg(input,option,keyword,err_string)
        ! the Q_max expression of Iyer et al. is posed in m^2/g
        internal_units = 'm^2/g'
        call InputReadAndConvertUnits(input, &
                                      new_pool%specific_surface_area, &
                                      internal_units,keyword,option)
      case('KD_LSCIE')
        call InputReadDouble(input,option,new_pool%Kd_lscie)
        call InputErrorMsg(input,option,keyword,err_string)
        ! the K_eq expression of Iyer et al. is posed in mL/g
        internal_units = 'mL/g'
        call InputReadAndConvertUnits(input,new_pool%Kd_lscie, &
                                      internal_units,keyword,option)
      case('MINERAL_DENSITY')
        call InputReadDouble(input,option,new_pool%mineral_density)
        call InputErrorMsg(input,option,keyword,err_string)
        ! paired with KD_LSCIE in K_eq, so likewise g/mL
        internal_units = 'g/mL'
        call InputReadAndConvertUnits(input,new_pool%mineral_density, &
                                      internal_units,keyword,option)
      case default
        call InputKeywordUnrecognized(input,keyword,err_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (len_trim(new_pool%species_name) == 0) then
    option%io_buffer = 'Q_SPECIES_NAME not defined in ' // trim(err_string)
    call PrintErrMsg(option)
  endif

  if (.not.associated(this%q_pool_list)) then
    this%q_pool_list => new_pool
  else
    cur_pool => this%q_pool_list
    do
      if (.not.associated(cur_pool%next)) exit
      cur_pool => cur_pool%next
    enddo
    cur_pool%next => new_pool
  endif
  nullify(new_pool)

end subroutine CarbonMENDnQReadSorptionSite

! ************************************************************************** !

subroutine CarbonMENDnQSetup(this,reaction,option)
  !
  ! Configures the reaction and associated data structures
  !
  ! Author: Richard Tran Mills
  ! Date: 10/26/26
  !
  use Material_Aux_module
  use Option_module
  use Reaction_Aux_module
  use Reaction_Immobile_Aux_module
  use String_module

  class(carbon_sandbox_mend_nq_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  type(carbon_sandbox_q_pool_type), pointer :: cur_pool
  PetscInt :: i
  PetscBool :: derived
  PetscReal :: Q_max_mg_per_g
  PetscReal :: sum_Q_max_mg_per_g
  PetscReal :: sum_mass_frac_over_density
  PetscReal :: sum_mass_fraction
  PetscReal, pointer :: Q_max_mg_per_g_array(:)

  call CarbonBaseSetup(this,reaction,option)

  this%B_species_index = &
    ReactionAuxGetPriSpecIDFromName(this%B_species_name,reaction,option)
  this%D_species_index = &
    ReactionAuxGetPriSpecIDFromName(this%D_species_name,reaction,option)
  this%EM_species_index = &
    ReactionAuxGetPriSpecIDFromName(this%EM_species_name,reaction,option)
  this%EP_species_index = &
    ReactionAuxGetPriSpecIDFromName(this%EP_species_name,reaction,option)
  this%IC_species_index = &
    ReactionAuxGetPriSpecIDFromName(this%IC_species_name,reaction,option)
  this%M_species_index = &
    ReactionImGetSpeciesIDFromName(this%M_species_name,reaction%immobile, &
                                   PETSC_FALSE,option) + &
    reaction%offset_immobile
  if (this%M_species_index <= 0) then
    option%io_buffer = 'Species "' // trim(this%M_species_name) // &
      '" not found among the available immobile species. Mineral associated &
      &carbon must be an immobile species.'
    call PrintErrMsg(option)
  endif
  this%P_species_index = &
    ReactionAuxGetPriSpecIDFromName(this%P_species_name,reaction,option)

  call CarbonMENDCheckCommonParams(this,option)
  call CarbonMENDSetOptionalDefaults(this,option)

  ! count the sorbed carbon pools
  this%n_q_pool = 0
  cur_pool => this%q_pool_list
  do
    if (.not.associated(cur_pool)) exit
    this%n_q_pool = this%n_q_pool + 1
    cur_pool => cur_pool%next
  enddo
  if (this%n_q_pool == 0) then
    option%io_buffer = 'At least one SORPTION_SITE must be defined in the &
      &MEND_NQ carbon sandbox.'
    call PrintErrMsg(option)
  endif

  allocate(this%Q_pool_index(this%n_q_pool))
  allocate(this%Q_pool_max(this%n_q_pool))
  allocate(this%Q_pool_kads(this%n_q_pool))
  allocate(this%Q_pool_kdes(this%n_q_pool))
  allocate(this%Q_pool_Keq_coef(this%n_q_pool))
  allocate(this%Q_pool_derived(this%n_q_pool))
  allocate(Q_max_mg_per_g_array(this%n_q_pool))
  this%Q_pool_index = UNINITIALIZED_INTEGER
  this%Q_pool_max = UNINITIALIZED_DOUBLE
  this%Q_pool_kads = UNINITIALIZED_DOUBLE
  this%Q_pool_kdes = UNINITIALIZED_DOUBLE
  this%Q_pool_Keq_coef = 0.d0
  this%Q_pool_derived = PETSC_FALSE
  Q_max_mg_per_g_array = 0.d0

  ! first pass: species indices, direct parameters, and Q_max for derived pools
  sum_Q_max_mg_per_g = 0.d0
  sum_mass_frac_over_density = 0.d0
  sum_mass_fraction = 0.d0
  i = 0
  cur_pool => this%q_pool_list
  do
    if (.not.associated(cur_pool)) exit
    i = i + 1

    this%Q_pool_index(i) = &
      ReactionImGetSpeciesIDFromName(cur_pool%species_name, &
                                     reaction%immobile,PETSC_FALSE,option) + &
      reaction%offset_immobile
    if (this%Q_pool_index(i) <= reaction%offset_immobile) then
      option%io_buffer = 'Species "' // trim(cur_pool%species_name) // &
        '" in SORPTION_SITE "' // trim(cur_pool%site_name) // '" was not &
        &found among the available immobile species. Each sorbed carbon pool &
        &must be an immobile species.'
      call PrintErrMsg(option)
    endif

    derived = .not.Uninitialized(cur_pool%mass_fraction)
    this%Q_pool_derived(i) = derived

    if (derived) then
      if (Initialized(cur_pool%Q_max) .or. &
          Initialized(cur_pool%K_ads) .or. &
          Initialized(cur_pool%K_des)) then
        option%io_buffer = 'SORPTION_SITE "' // trim(cur_pool%site_name) // &
          '" mixes mineral properties with Q_MAX/K_ADS/K_DES. A sorption &
          &site must specify either the mineral properties or the sorption &
          &parameters, not both.'
        call PrintErrMsg(option)
      endif
      call CarbonMENDnQCheckMineral(cur_pool,this,option)
      ! Q_max,i = w*SA*SD*MW_OM*f_C*1.d21/N_A, in mg C/g soil
      Q_max_mg_per_g = cur_pool%mass_fraction * &
                       cur_pool%specific_surface_area * &
                       cur_pool%site_density * this%MW_OM * this%f_C * &
                       1.d21 / AVOGADRO_NUMBER
      Q_max_mg_per_g_array(i) = Q_max_mg_per_g
      sum_Q_max_mg_per_g = sum_Q_max_mg_per_g + Q_max_mg_per_g
      sum_mass_frac_over_density = sum_mass_frac_over_density + &
                                   cur_pool%mass_fraction / &
                                   cur_pool%mineral_density
      sum_mass_fraction = sum_mass_fraction + cur_pool%mass_fraction
      this%Q_pool_max(i) = Q_max_mg_per_g
    else
      if (Uninitialized(cur_pool%Q_max)) then
        option%io_buffer = 'Q_MAX not defined for SORPTION_SITE "' // &
          trim(cur_pool%site_name) // '" in the MEND_NQ carbon sandbox.'
        call PrintErrMsg(option)
      endif
      if (Uninitialized(cur_pool%K_ads)) then
        option%io_buffer = 'K_ADS not defined for SORPTION_SITE "' // &
          trim(cur_pool%site_name) // '" in the MEND_NQ carbon sandbox.'
        call PrintErrMsg(option)
      endif
      if (Uninitialized(cur_pool%K_des)) then
        option%io_buffer = 'K_DES not defined for SORPTION_SITE "' // &
          trim(cur_pool%site_name) // '" in the MEND_NQ carbon sandbox.'
        call PrintErrMsg(option)
      endif
      this%Q_pool_max(i) = cur_pool%Q_max
      this%Q_pool_kads(i) = cur_pool%K_ads
      this%Q_pool_kdes(i) = cur_pool%K_des
    endif
    cur_pool => cur_pool%next
  enddo

  ! second pass: the derived rate constants, which need the sums above
  if (sum_Q_max_mg_per_g > 0.d0) then
    if (Uninitialized(this%K_ads_base)) then
      option%io_buffer = 'K_ADS_BASE must be defined in the MINERALOGY block &
        &of the MEND_NQ carbon sandbox when any SORPTION_SITE is defined by &
        &mineral properties.'
      call PrintErrMsg(option)
    endif
    if (dabs(sum_mass_fraction-1.d0) > 1.d-6) then
      option%io_buffer = 'The mass fractions of the mineral-derived &
        &SORPTION_SITEs in the MEND_NQ carbon sandbox sum to ' // &
        trim(StringWrite(sum_mass_fraction)) // ' rather than 1.0.'
      call PrintWrnMsg(option)
    endif
    i = 0
    cur_pool => this%q_pool_list
    do
      if (.not.associated(cur_pool)) exit
      i = i + 1
      if (this%Q_pool_derived(i)) then
        ! k_ads,i partitions K_ads_base by each mineral's share of the
        ! total sorption capacity, which is its share of reactive surface
        this%Q_pool_kads(i) = this%K_ads_base * &
                              Q_max_mg_per_g_array(i) / sum_Q_max_mg_per_g
        ! K_eq,i less the (1-phi)/(phi*S_w) factor, which is applied in
        ! Evaluate because porosity and saturation are live state
        this%Q_pool_Keq_coef(i) = cur_pool%Kd_lscie * &
                                  cur_pool%mass_fraction / &
                                  (sum_mass_frac_over_density * &
                                   Q_max_mg_per_g_array(i))
      endif
      cur_pool => cur_pool%next
    enddo
  endif

  deallocate(Q_max_mg_per_g_array)
  nullify(Q_max_mg_per_g_array)

  call this%EnforceConcentrationUnits(CARBON_UNITS_MOLE_PER_KG_SOIL,option)

end subroutine CarbonMENDnQSetup

! ************************************************************************** !

subroutine CarbonMENDnQCheckMineral(q_pool,this,option)
  !
  ! Verifies that a mineral-derived sorption site is fully specified
  !
  ! Author: Richard Tran Mills
  ! Date: 10/26/26
  !
  use Option_module

  type(carbon_sandbox_q_pool_type) :: q_pool
  class(carbon_sandbox_mend_nq_type) :: this
  type(option_type) :: option

  call CheckMineralParam(q_pool%specific_surface_area, &
                         'SPECIFIC_SURFACE_AREA')
  call CheckMineralParam(q_pool%site_density,'SITE_DENSITY')
  call CheckMineralParam(q_pool%Kd_lscie,'KD_LSCIE')
  call CheckMineralParam(q_pool%mineral_density,'MINERAL_DENSITY')
  call CheckMineralogyParam(this%MW_OM,'MW_OM')
  call CheckMineralogyParam(this%f_C,'F_C')

contains

subroutine CheckMineralParam(var,varname)

  PetscReal :: var
  character(len=*) :: varname

  if (Uninitialized(var)) then
    option%io_buffer = trim(varname) // ' not defined for SORPTION_SITE "' // &
      trim(q_pool%site_name) // '" in the MEND_NQ carbon sandbox.'
    call PrintErrMsg(option)
  endif

end subroutine CheckMineralParam

subroutine CheckMineralogyParam(var,varname)

  PetscReal :: var
  character(len=*) :: varname

  if (Uninitialized(var)) then
    option%io_buffer = trim(varname) // ' not defined in the MINERALOGY &
      &block of the MEND_NQ carbon sandbox, but SORPTION_SITE "' // &
      trim(q_pool%site_name) // '" is defined by mineral properties.'
    call PrintErrMsg(option)
  endif

end subroutine CheckMineralogyParam

end subroutine CarbonMENDnQCheckMineral

! ************************************************************************** !

subroutine CarbonMENDnQEvaluate(this,Residual,Jacobian,compute_derivative, &
                                rt_auxvar,global_auxvar,material_auxvar, &
                                reaction,option)
  !
  ! Evaluates the rate expression
  !
  ! Author: Richard Tran Mills
  ! Date: 10/26/26
  !
  use Global_Aux_module
  use Material_Aux_module
  use Option_module
  use Reaction_Aux_module
  use Reaction_Inhibition_Aux_module
  use Reactive_Transport_Aux_module
  use Utility_module

  class(carbon_sandbox_mend_nq_type) :: this
  class(reaction_rt_type) :: reaction
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(:)
  PetscReal :: Jacobian(:,:)
  PetscBool :: compute_derivative
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option

  PetscReal :: B, D, EM, EP, IC, M, P, Q
  PetscReal :: D_monod, M_monod, P_monod
  PetscReal :: one_over_E_C
  PetscReal :: F1, F2, F3, F4, F5, F6, F7, F8, F9EP, F9EM, F10EP, F10EM
  PetscReal :: dB, dD, dEM, dEP, dIC, dM, dP, dQ
  PetscReal :: conc_units_conversion
  PetscReal :: V_P_adj, V_M_adj, V_D_adj, K_P_adj, K_M_adj, K_D_adj
  PetscReal :: m_R_adj, E_C_adj
  PetscReal :: K_ads_adj, K_des_adj
  PetscReal :: K_ads_scale, K_des_scale
  PetscReal :: K_des_i, K_eq_i, saturation_factor
  PetscReal :: t
  PetscReal :: tempreal
  PetscInt :: i

  call this%MapStateVariables(rt_auxvar,global_auxvar,material_auxvar, &
                              reaction,option)

  t = this%aux%temperature
  ! mgC/gsoil = molC/kgsoil*gC/molC*mgC/gC*kgsoil/gsoil
  conc_units_conversion = 30.d0*1000.d0*1.d-3
  B = this%aux%conc(this%B_species_index)*conc_units_conversion
  D = this%aux%conc(this%D_species_index)*conc_units_conversion
  EM = this%aux%conc(this%EM_species_index)*conc_units_conversion
  EP = this%aux%conc(this%EP_species_index)*conc_units_conversion
  IC = this%aux%conc(this%IC_species_index)*conc_units_conversion
  M = this%aux%conc(this%M_species_index)*conc_units_conversion
  P = this%aux%conc(this%P_species_index)*conc_units_conversion

  ! activation energies are from Wang et al. (2013), Table 3, scaled from the
  ! kJ/mol reported there to the J/mol expected by Arrhenius()
  V_P_adj = this%V_P*Arrhenius(53.d3,t,this%reference_temperature)
  tempreal = Arrhenius(47.d3,t,this%reference_temperature)
  V_M_adj = this%V_M*tempreal
  V_D_adj = this%V_D*tempreal
  tempreal = Arrhenius(30.d3,t,this%reference_temperature)
  K_P_adj = this%K_P*tempreal
  K_M_adj = this%K_M*tempreal
  K_D_adj = this%K_D*tempreal
  m_R_adj = this%m_R*Arrhenius(20.d3,t,this%reference_temperature)
  K_ads_scale = Arrhenius(5.d3,t,this%reference_temperature)
  K_des_scale = Arrhenius(20.d3,t,this%reference_temperature)
  E_C_adj = &
    max(min(this%E_C - 0.012d0*(t-this%reference_temperature),1.d0),0.d0)

  D_monod = D/(K_D_adj+D)
  M_monod = M/(K_M_adj+M)
  P_monod = P/(K_P_adj+P)
  one_over_E_C = 1.d0 / E_C_adj
  F1 = one_over_E_C*(V_D_adj+m_R_adj)*B*D_monod
  F2 = V_P_adj*EP*P_monod
  F3 = V_M_adj*EM*M_monod
  F4 = (one_over_E_C-1.d0)*V_D_adj*B*D_monod
  F5 = (one_over_E_C-1.d0)*m_R_adj*B*D_monod
  F8 = (1.d0-this%p_EP-this%p_EM)*m_R_adj*B
  F9EM = this%p_EM*m_R_adj*B
  F9EP = this%p_EP*m_R_adj*B
  F10EM = this%r_EM*EM
  F10EP = this%r_EP*EP

  ! (1-phi)/(phi*S_w) in K_eq, from the live flow state. Guard against a dry
  ! cell, where desorption is switched off rather than allowed to blow up.
  saturation_factor = this%aux%porosity*this%aux%liquid_saturation
  if (saturation_factor > 0.d0) then
    saturation_factor = (1.d0-this%aux%porosity)/saturation_factor
  else
    saturation_factor = 0.d0
  endif

  ! sum the adsorption and desorption fluxes over the sorbed carbon pools,
  ! writing each pool's residual as we go
  F6 = 0.d0
  F7 = 0.d0
  do i = 1, this%n_q_pool
    Q = this%aux%conc(this%Q_pool_index(i))*conc_units_conversion
    K_ads_adj = this%Q_pool_kads(i)*K_ads_scale
    if (this%Q_pool_derived(i)) then
      K_eq_i = this%Q_pool_Keq_coef(i)*saturation_factor
      if (K_eq_i > 0.d0) then
        K_des_i = this%Q_pool_kads(i)/K_eq_i
      else
        K_des_i = 0.d0
      endif
    else
      K_des_i = this%Q_pool_kdes(i)
    endif
    K_des_adj = K_des_i*K_des_scale
    dQ = K_ads_adj*(1.d0-Q/this%Q_pool_max(i))*D - &
         K_des_adj*Q/this%Q_pool_max(i)
    F6 = F6 + K_ads_adj*(1.d0-Q/this%Q_pool_max(i))*D
    F7 = F7 + K_des_adj*Q/this%Q_pool_max(i)
    Residual(this%Q_pool_index(i)) = Residual(this%Q_pool_index(i)) - dQ
  enddo

  dB = F1 - (F4 + F5) - F8 - (F9EP + F9EM)
  dD = this%I_D + this%f_D*F2 + this%g_D*F8 + F3 + (F10EP + F10EM) - &
       F1 - (F6 - F7)
  dEM = F9EM - F10EM
  dEP = F9EP - F10EP
  dM = (1.d0-this%f_D)*F2 - F3
  dP = this%I_P + (1.d0-this%g_D)*F8 - F2
  dIC = F4 + F5

  Residual(this%B_species_index) = Residual(this%B_species_index) - dB
  Residual(this%D_species_index) = Residual(this%D_species_index) - dD
  Residual(this%EM_species_index) = Residual(this%EM_species_index) - dEM
  Residual(this%EP_species_index) = Residual(this%EP_species_index) - dEP
  Residual(this%IC_species_index) = Residual(this%IC_species_index) - dIC
  Residual(this%M_species_index) = Residual(this%M_species_index) - dM
  Residual(this%P_species_index) = Residual(this%P_species_index) - dP

end subroutine CarbonMENDnQEvaluate

! ************************************************************************** !

subroutine CarbonMENDnQStrip(this)
  !
  ! Destroys members of the carbon sandbox MEND_NQ object
  !
  ! Author: Richard Tran Mills
  ! Date: 10/26/26
  !
  use Utility_module

  class(carbon_sandbox_mend_nq_type) :: this

  type(carbon_sandbox_q_pool_type), pointer :: cur_pool, next_pool

  call CarbonBaseStrip(this)

  cur_pool => this%q_pool_list
  do
    if (.not.associated(cur_pool)) exit
    next_pool => cur_pool%next
    deallocate(cur_pool)
    cur_pool => next_pool
  enddo
  nullify(this%q_pool_list)

  call DeallocateArray(this%Q_pool_index)
  call DeallocateArray(this%Q_pool_max)
  call DeallocateArray(this%Q_pool_kads)
  call DeallocateArray(this%Q_pool_kdes)
  call DeallocateArray(this%Q_pool_Keq_coef)
  if (associated(this%Q_pool_derived)) then
    deallocate(this%Q_pool_derived)
    nullify(this%Q_pool_derived)
  endif

end subroutine CarbonMENDnQStrip

end module Carbon_Sandbox_MEND_nQ_class
