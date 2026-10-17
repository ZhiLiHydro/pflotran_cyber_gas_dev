module Adjoint_Param_module

#include "petsc/finclude/petscvec.h"
  use petscvec
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

  use PFLOTRAN_Constants_module

  ! The ADJOINT block's parameters, in deck order. A parameter is data: its
  ! kind, the reaction, Monod or inhibition term, mineral or material it
  ! names, and its value; the adjoint of the equation model differentiates
  ! the model with respect to it (Adjoint_RT_DAE_class).
  !
  !   PARAMETER RATE_CONSTANT            KINETICS MICROBIAL, REACTION_INDEX
  !                                      or KINETICS MINERAL, MINERAL_NAME
  !   PARAMETER HALF_SATURATION_CONSTANT KINETICS MICROBIAL, REACTION_INDEX,
  !                                      SPECIES_NAME (the Monod term)
  !   PARAMETER INHIBITION_CONSTANT      same, the inhibition term (MONOD or
  !                                      SMOOTHSTEP): its threshold magnitude
  !   PARAMETER SMOOTHSTEP_INTERVAL      same, a SMOOTHSTEP term's log10
  !                                      interval
  !   PARAMETER YIELD                    KINETICS MICROBIAL, REACTION_INDEX
  !                                      (a reaction with a BIOMASS block)
  !   PARAMETER ACTIVATION_ENERGY        KINETICS MICROBIAL, REACTION_INDEX
  !                                      (a reaction with a positive
  !                                      ACTIVATION_ENERGY) [J/mol]
  !   PARAMETER DECAY_RATE_CONSTANT      REACTION_INDEX (of the
  !                                      IMMOBILE_DECAY_REACTION blocks)
  !   PARAMETER PERMEABILITY             MATERIAL (isotropic; scales the
  !                                      material's cell permeabilities)
  !   PARAMETER INITIAL_CONCENTRATION    SPECIES_NAME, REGION: the initial
  !                                      free molality [mol/kg] of a primary
  !                                      species, or concentration
  !                                      [mol/m^3 bulk] of an immobile one,
  !                                      in the region's cells (one value);
  !                                      not a model parameter: the adjoint
  !                                      sets its initial state
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26

  implicit none

  private

  PetscInt, parameter, public :: ADJ_PARAM_RATE_CONSTANT = 1
  PetscInt, parameter, public :: ADJ_PARAM_HALF_SATURATION = 2
  PetscInt, parameter, public :: ADJ_PARAM_INHIBITION_CONSTANT = 3
  PetscInt, parameter, public :: ADJ_PARAM_SMOOTHSTEP_INTERVAL = 4
  PetscInt, parameter, public :: ADJ_PARAM_YIELD = 5
  PetscInt, parameter, public :: ADJ_PARAM_DECAY_RATE = 6
  PetscInt, parameter, public :: ADJ_PARAM_PERMEABILITY = 7
  PetscInt, parameter, public :: ADJ_PARAM_INITIAL_CONCENTRATION = 8
  PetscInt, parameter, public :: ADJ_PARAM_ACTIVATION_ENERGY = 9

  PetscInt, parameter, public :: ADJ_KINETICS_MICROBIAL = 1
  PetscInt, parameter, public :: ADJ_KINETICS_MINERAL = 2

  type, public :: adjoint_param_type
    PetscInt :: id = UNINITIALIZED_INTEGER      ! position in list
    PetscInt :: kind = UNINITIALIZED_INTEGER
    character(len=MAXSTRINGLENGTH) :: name = ''   ! label for output
    PetscReal :: value = UNINITIALIZED_DOUBLE
    PetscInt :: kinetics = UNINITIALIZED_INTEGER
    ! microbial or immobile-decay reaction (1-based, in deck order)
    PetscInt :: reaction_index = UNINITIALIZED_INTEGER
    ! the species selecting a Monod or inhibition term, and the term
    ! (resolved at Setup)
    character(len=MAXWORDLENGTH) :: species_name = ''
    PetscInt :: term = UNINITIALIZED_INTEGER
    ! INHIBITION_CONSTANT: sign of inhibition_C (negative: inhibit above)
    PetscReal :: threshold_sign = 1.d0
    ! RATE_CONSTANT KINETICS MINERAL
    character(len=MAXWORDLENGTH) :: mineral_name = ''
    PetscInt :: imnrl = 0
    ! PERMEABILITY: the material, its value and cell permeabilities at Setup
    character(len=MAXWORDLENGTH) :: material_name = ''
    PetscInt :: material_id = UNINITIALIZED_INTEGER
    PetscReal :: value0 = 0.d0
    PetscReal, allocatable :: perm0(:,:)
    ! INITIAL_CONCENTRATION: the region, the species' entry in a cell's c
    ! coordinates, and the region's active local cells (at Setup)
    character(len=MAXWORDLENGTH) :: region_name = ''
    PetscInt :: component = 0
    PetscInt, allocatable :: cells(:)
    type(adjoint_param_type), pointer :: next => null()
  contains
    procedure, public :: Setup => AdjParamSetup
    procedure, public :: SetValue => AdjParamSetValue
  end type adjoint_param_type

  type(adjoint_param_type), pointer, public :: adj_param_list => null()

  public :: AdjParamInit, &
            AdjParamRead, &
            AdjParamCount, &
            AdjParamDestroyList

contains

! ************************************************************************** !

subroutine AdjParamInit()
  !
  ! Initializes the parameter list, destroying any existing one.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  implicit none

  if (associated(adj_param_list)) then
    call AdjParamDestroyList()
  endif
  nullify(adj_param_list)

end subroutine AdjParamInit

! ************************************************************************** !

subroutine AdjParamRead(input,option)
  !
  ! Reads one PARAMETER <class> block and appends the parameter to the
  ! list.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  use Option_module
  use String_module
  use Input_Aux_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: class_keyword, keyword, word
  character(len=MAXSTRINGLENGTH) :: error_string, label
  type(adjoint_param_type), pointer :: new_param, cur_param
  PetscBool :: allowed

  error_string = 'SUBSURFACE,ADJOINT,PARAMETER'

  call InputReadCard(input,option,class_keyword)
  call InputErrorMsg(input,option,'class keyword',error_string)
  call StringToUpper(class_keyword)

  allocate(new_param)
  select case(trim(class_keyword))
    case('RATE_CONSTANT')
      new_param%kind = ADJ_PARAM_RATE_CONSTANT
    case('HALF_SATURATION_CONSTANT')
      new_param%kind = ADJ_PARAM_HALF_SATURATION
    case('INHIBITION_CONSTANT')
      new_param%kind = ADJ_PARAM_INHIBITION_CONSTANT
    case('SMOOTHSTEP_INTERVAL')
      new_param%kind = ADJ_PARAM_SMOOTHSTEP_INTERVAL
    case('YIELD')
      new_param%kind = ADJ_PARAM_YIELD
    case('DECAY_RATE_CONSTANT')
      new_param%kind = ADJ_PARAM_DECAY_RATE
    case('PERMEABILITY')
      new_param%kind = ADJ_PARAM_PERMEABILITY
    case('INITIAL_CONCENTRATION')
      new_param%kind = ADJ_PARAM_INITIAL_CONCENTRATION
    case('ACTIVATION_ENERGY')
      new_param%kind = ADJ_PARAM_ACTIVATION_ENERGY
    case default
      deallocate(new_param)
      call InputKeywordUnrecognized(input,class_keyword,error_string,option)
      return
  end select
  error_string = trim(error_string) // ',' // trim(class_keyword)

  input%ierr = INPUT_ERROR_NONE
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    ! the keywords of each class
    select case(trim(keyword))
      case('KINETICS')
        allowed = new_param%kind /= ADJ_PARAM_DECAY_RATE .and. &
                  new_param%kind /= ADJ_PARAM_PERMEABILITY .and. &
                  new_param%kind /= ADJ_PARAM_INITIAL_CONCENTRATION
      case('REACTION_INDEX')
        allowed = new_param%kind /= ADJ_PARAM_PERMEABILITY .and. &
                  new_param%kind /= ADJ_PARAM_INITIAL_CONCENTRATION
      case('SPECIES_NAME')
        allowed = new_param%kind == ADJ_PARAM_HALF_SATURATION .or. &
                  new_param%kind == ADJ_PARAM_INHIBITION_CONSTANT .or. &
                  new_param%kind == ADJ_PARAM_SMOOTHSTEP_INTERVAL .or. &
                  new_param%kind == ADJ_PARAM_INITIAL_CONCENTRATION
      case('MINERAL_NAME')
        allowed = new_param%kind == ADJ_PARAM_RATE_CONSTANT
      case('MATERIAL')
        allowed = new_param%kind == ADJ_PARAM_PERMEABILITY
      case('REGION')
        allowed = new_param%kind == ADJ_PARAM_INITIAL_CONCENTRATION
      case default
        allowed = PETSC_FALSE
    end select
    if (.not.allowed) &
      call InputKeywordUnrecognized(input,keyword,error_string,option)

    select case(trim(keyword))
      case('KINETICS')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)
        call StringToUpper(word)
        select case(trim(word))
          case('MICROBIAL')
            new_param%kinetics = ADJ_KINETICS_MICROBIAL
          case('MINERAL')
            select case(new_param%kind)
              case(ADJ_PARAM_RATE_CONSTANT)
                new_param%kinetics = ADJ_KINETICS_MINERAL
              case(ADJ_PARAM_YIELD)
                option%io_buffer = 'KINETICS MINERAL is not supported for &
                  &the YIELD adjoint parameter (yield is a microbial-biomass &
                  &concept). Use KINETICS MICROBIAL.'
                call PrintErrMsg(option)
              case default
                option%io_buffer = 'KINETICS must be MICROBIAL in ' // &
                                   trim(error_string) // '.'
                call PrintErrMsg(option)
            end select
          case default
            if (new_param%kind == ADJ_PARAM_RATE_CONSTANT) then
              option%io_buffer = 'KINETICS must be MICROBIAL or MINERAL in ' &
                                 // trim(error_string) // '.'
            else
              option%io_buffer = 'KINETICS must be MICROBIAL in ' // &
                                 trim(error_string) // '.'
            endif
            call PrintErrMsg(option)
        end select
      case('REACTION_INDEX')
        call InputReadInt(input,option,new_param%reaction_index)
        call InputErrorMsg(input,option,keyword,error_string)
      case('SPECIES_NAME')
        call InputReadWord(input,option,new_param%species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)
      case('MINERAL_NAME')
        call InputReadWord(input,option,new_param%mineral_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)
      case('MATERIAL')
        call InputReadWord(input,option,new_param%material_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)
      case('REGION')
        call InputReadWord(input,option,new_param%region_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)
    end select
  enddo
  call InputPopBlock(input,option)

  select case(new_param%kind)
    case(ADJ_PARAM_RATE_CONSTANT)
      select case(new_param%kinetics)
        case(ADJ_KINETICS_MICROBIAL)
          if (Uninitialized(new_param%reaction_index)) then
            option%io_buffer = 'REACTION_INDEX must be specified for &
              &RATE_CONSTANT KINETICS MICROBIAL.'
            call PrintErrMsg(option)
          endif
          if (len_trim(new_param%mineral_name) > 0) then
            option%io_buffer = 'MINERAL_NAME does not apply to RATE_CONSTANT &
              &KINETICS MICROBIAL.'
            call PrintErrMsg(option)
          endif
          write(new_param%name,'(a,"[",i0,"]")') &
                'MICROBIAL_RATE_CONSTANT', new_param%reaction_index
        case(ADJ_KINETICS_MINERAL)
          if (len_trim(new_param%mineral_name) == 0) then
            option%io_buffer = 'MINERAL_NAME must be specified for &
              &RATE_CONSTANT KINETICS MINERAL.'
            call PrintErrMsg(option)
          endif
          if (Initialized(new_param%reaction_index)) then
            option%io_buffer = 'REACTION_INDEX does not apply to &
              &RATE_CONSTANT KINETICS MINERAL; name the mineral with &
              &MINERAL_NAME.'
            call PrintErrMsg(option)
          endif
          new_param%name = 'MINERAL_RATE_CONSTANT[' // &
                           trim(new_param%mineral_name) // ']'
        case default
          option%io_buffer = 'KINETICS must be specified for RATE_CONSTANT &
            &parameter.'
          call PrintErrMsg(option)
      end select
    case(ADJ_PARAM_HALF_SATURATION,ADJ_PARAM_INHIBITION_CONSTANT, &
         ADJ_PARAM_SMOOTHSTEP_INTERVAL)
      if (Uninitialized(new_param%kinetics)) then
        option%io_buffer = 'KINETICS must be specified for ' // &
          trim(class_keyword) // ' parameter.'
        call PrintErrMsg(option)
      endif
      if (Uninitialized(new_param%reaction_index)) then
        option%io_buffer = 'REACTION_INDEX must be specified for ' // &
          trim(class_keyword) // ' parameter.'
        call PrintErrMsg(option)
      endif
      if (len_trim(new_param%species_name) == 0) then
        option%io_buffer = 'SPECIES_NAME must be specified for ' // &
          trim(class_keyword) // ' parameter.'
        call PrintErrMsg(option)
      endif
      select case(new_param%kind)
        case(ADJ_PARAM_HALF_SATURATION)
          word = 'HALF_SAT'
        case(ADJ_PARAM_INHIBITION_CONSTANT)
          word = 'INHIB_C'
        case(ADJ_PARAM_SMOOTHSTEP_INTERVAL)
          word = 'SS_INT'
      end select
      write(label,'(a,"[",i0,",",a,"]")') &
            trim(word), new_param%reaction_index, trim(new_param%species_name)
      new_param%name = label(1:len(new_param%name))
    case(ADJ_PARAM_YIELD,ADJ_PARAM_ACTIVATION_ENERGY)
      if (Uninitialized(new_param%kinetics)) then
        option%io_buffer = 'KINETICS must be specified for ' // &
          trim(class_keyword) // ' parameter.'
        call PrintErrMsg(option)
      endif
      if (Uninitialized(new_param%reaction_index)) then
        option%io_buffer = 'REACTION_INDEX must be specified for ' // &
          trim(class_keyword) // ' KINETICS MICROBIAL.'
        call PrintErrMsg(option)
      endif
      write(new_param%name,'(a,"_",a,"[",i0,"]")') 'MICROBIAL', &
            trim(class_keyword), new_param%reaction_index
    case(ADJ_PARAM_DECAY_RATE)
      if (Uninitialized(new_param%reaction_index)) then
        option%io_buffer = 'REACTION_INDEX must be specified for &
          &DECAY_RATE_CONSTANT.'
        call PrintErrMsg(option)
      endif
      write(new_param%name,'(a,"[",i0,"]")') &
            'IMMOBILE_DECAY_RATE', new_param%reaction_index
    case(ADJ_PARAM_PERMEABILITY)
      if (len_trim(new_param%material_name) == 0) then
        call PrintErrMsg(option,'PARAMETER PERMEABILITY requires MATERIAL.')
      endif
      new_param%name = 'PERMEABILITY[' // trim(new_param%material_name) // ']'
    case(ADJ_PARAM_INITIAL_CONCENTRATION)
      if (len_trim(new_param%species_name) == 0 .or. &
          len_trim(new_param%region_name) == 0) then
        call PrintErrMsg(option,'PARAMETER INITIAL_CONCENTRATION requires &
          &SPECIES_NAME and REGION.')
      endif
      new_param%name = 'INITIAL_CONCENTRATION[' // &
        trim(new_param%species_name) // ',' // &
        trim(new_param%region_name) // ']'
  end select

  if (.not.associated(adj_param_list)) then
    adj_param_list => new_param
    new_param%id = 1
  else
    ! the name identifies the parameter: a second copy would be a second
    ! optimization variable that only one value of reaches the model
    cur_param => adj_param_list
    do
      if (cur_param%name == new_param%name) then
        call PrintErrMsg(option,'SUBSURFACE,ADJOINT: PARAMETER ' // &
          trim(new_param%name) // ' is given twice.')
      endif
      if (.not.associated(cur_param%next)) exit
      cur_param => cur_param%next
    enddo
    cur_param%next => new_param
    new_param%id = cur_param%id + 1
  endif

end subroutine AdjParamRead

! ************************************************************************** !

subroutine AdjParamSetup(this,realization,option)
  !
  ! Resolves the reaction, term, mineral or material and reads the
  ! parameter's value.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  use Option_module
  use Realization_Subsurface_class
  use Reaction_Mineral_Aux_module, only : mineral_type, &
                                          ReactionMnrlGetMnrlIDFromName
  use Reaction_Microbial_Aux_module, only : microbial_type
  use Reaction_Inhibition_Aux_module, only : INHIBITION_MONOD, &
                                             INHIBITION_SMOOTHSTEP
  use Material_module, only : material_property_type, &
                              MaterialPropGetPtrFromList
  use Patch_module
  use Field_module
  use Region_module, only : region_type, RegionGetPtrFromList
  use Reaction_Aux_module, only : reaction_rt_type, &
                                  ReactionAuxGetPriSpecIDFromName
  use Reaction_Immobile_Aux_module, only : ReactionImGetSpeciesIDFromName

  implicit none

  class(adjoint_param_type) :: this
  class(realization_subsurface_type) :: realization
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: class_keyword
  type(mineral_type), pointer :: mineral
  type(microbial_type), pointer :: microbial
  type(material_property_type), pointer :: mat
  type(patch_type), pointer :: patch
  type(region_type), pointer :: region
  class(reaction_rt_type), pointer :: reaction
  PetscReal, pointer :: xx_p(:)
  PetscReal :: interval, threshold, extremes(2)
  PetscInt, allocatable :: cells(:)
  PetscInt :: local_id, ncells, i, id
  PetscErrorCode :: ierr

  select case(this%kind)
    case(ADJ_PARAM_RATE_CONSTANT)
      select case(this%kinetics)
        case(ADJ_KINETICS_MICROBIAL)
          if (realization%reaction%microbial%nrxn == 0) then
            option%io_buffer = 'RATE_CONSTANT KINETICS MICROBIAL specified, &
              &but the chemistry deck contains no MICROBIAL_REACTION block.'
            call PrintErrMsg(option)
          endif
          if (this%reaction_index < 1 .or. &
              this%reaction_index > realization%reaction%microbial%nrxn) then
            option%io_buffer = 'REACTION_INDEX out of range for &
              &RATE_CONSTANT KINETICS MICROBIAL.'
            call PrintErrMsg(option)
          endif
          this%value = realization%reaction%microbial% &
                         rate_constant(this%reaction_index)
        case(ADJ_KINETICS_MINERAL)
          mineral => realization%reaction%mineral
          this%imnrl = ReactionMnrlGetMnrlIDFromName(this%mineral_name, &
                         mineral,PETSC_TRUE,PETSC_FALSE,option)
          if (this%imnrl <= 0) then
            option%io_buffer = 'RATE_CONSTANT KINETICS MINERAL: "' // &
              trim(this%mineral_name) // '" is not a kinetic mineral.'
            call PrintErrMsg(option)
          endif
          if (mineral%kinmnrl_precip_rate_constant(this%imnrl) /= &
              mineral%kinmnrl_dissol_rate_constant(this%imnrl)) then
            option%io_buffer = 'RATE_CONSTANT KINETICS MINERAL "' // &
              trim(this%mineral_name) // '" requires equal precipitation &
              &and dissolution rate constants (a single RATE_CONSTANT).'
            call PrintErrMsg(option)
          endif
          this%value = mineral%kinmnrl_dissol_rate_constant(this%imnrl)
      end select
      if (.not.ieee_is_finite(this%value) .or. this%value < 0.d0) then
        call PrintErrMsg(option,'RATE_CONSTANT adjoint requires a finite, &
          &nonnegative rate constant.')
      endif
    case(ADJ_PARAM_HALF_SATURATION,ADJ_PARAM_INHIBITION_CONSTANT, &
         ADJ_PARAM_SMOOTHSTEP_INTERVAL)
      select case(this%kind)
        case(ADJ_PARAM_HALF_SATURATION)
          class_keyword = 'HALF_SATURATION_CONSTANT'
        case(ADJ_PARAM_INHIBITION_CONSTANT)
          class_keyword = 'INHIBITION_CONSTANT'
        case(ADJ_PARAM_SMOOTHSTEP_INTERVAL)
          class_keyword = 'SMOOTHSTEP_INTERVAL'
      end select
      if (realization%reaction%microbial%nrxn == 0) then
        option%io_buffer = trim(class_keyword) // ' KINETICS MICROBIAL &
          &specified, but the chemistry deck contains no MICROBIAL_REACTION &
          &block.'
        call PrintErrMsg(option)
      endif
      microbial => realization%reaction%microbial
      if (this%reaction_index < 1 .or. &
          this%reaction_index > microbial%nrxn) then
        option%io_buffer = 'REACTION_INDEX out of range for ' // &
          trim(class_keyword) // ' KINETICS MICROBIAL.'
        call PrintErrMsg(option)
      endif
      call AdjParamSelectTerm(this,realization,option)
      select case(this%kind)
        case(ADJ_PARAM_HALF_SATURATION)
          this%value = microbial%monod_K(this%term)
          if (.not.ieee_is_finite(this%value) .or. this%value <= 0.d0) then
            call PrintErrMsg(option,'HALF_SATURATION_CONSTANT adjoint &
              &requires a positive half-saturation constant.')
          endif
        case(ADJ_PARAM_INHIBITION_CONSTANT)
          this%threshold_sign = sign(1.d0,microbial%inhibition_C(this%term))
          this%value = abs(microbial%inhibition_C(this%term))
          select case(microbial%inhibition_type(this%term))
            case(INHIBITION_MONOD,INHIBITION_SMOOTHSTEP)
            case default
              option%io_buffer = 'INHIBITION_CONSTANT (threshold derivative) &
                &supports TYPE MONOD and SMOOTHSTEP only; THRESHOLD (arctan) &
                &is not implemented here.'
              call PrintErrMsg(option)
          end select
          if (.not.ieee_is_finite(this%value) .or. this%value <= 0.d0) then
            call PrintErrMsg(option,'INHIBITION_CONSTANT adjoint requires a &
              &finite, positive threshold magnitude.')
          endif
          if (microbial%inhibition_type(this%term) == &
              INHIBITION_SMOOTHSTEP) then
            interval = microbial%inhibition_C2(this%term)
            if (.not.ieee_is_finite(interval) .or. interval <= 0.d0) then
              call PrintErrMsg(option,'SMOOTHSTEP adjoint requires a finite, &
                &positive interval.')
            endif
          endif
        case(ADJ_PARAM_SMOOTHSTEP_INTERVAL)
          threshold = microbial%inhibition_C(this%term)
          this%value = microbial%inhibition_C2(this%term)
          if (.not.ieee_is_finite(this%value) .or. this%value <= 0.d0) then
            call PrintErrMsg(option,'SMOOTHSTEP adjoint requires a finite, &
              &positive interval.')
          endif
          if (.not.ieee_is_finite(threshold) .or. abs(threshold) == 0.d0) then
            call PrintErrMsg(option,'SMOOTHSTEP adjoint requires a finite, &
              &nonzero threshold.')
          endif
      end select
    case(ADJ_PARAM_YIELD)
      if (realization%reaction%microbial%nrxn == 0) then
        option%io_buffer = 'YIELD KINETICS MICROBIAL specified, but the &
          &chemistry deck contains no MICROBIAL_REACTION block.'
        call PrintErrMsg(option)
      endif
      if (this%reaction_index < 1 .or. &
          this%reaction_index > realization%reaction%microbial%nrxn) then
        option%io_buffer = 'REACTION_INDEX out of range for &
          &YIELD KINETICS MICROBIAL.'
        call PrintErrMsg(option)
      endif
      ! biomassid == 0: no BIOMASS block, the yield is unused
      if (realization%reaction%microbial%biomassid(this%reaction_index) &
          == 0) then
        option%io_buffer = 'YIELD parameter requires the target &
          &MICROBIAL_REACTION to have a BIOMASS block, but REACTION_INDEX &
          &has none.'
        call PrintErrMsg(option)
      endif
      this%value = &
        realization%reaction%microbial%biomass_yield(this%reaction_index)
      if (.not.ieee_is_finite(this%value) .or. this%value < 0.d0) then
        call PrintErrMsg(option,'YIELD adjoint requires a finite, &
          &nonnegative parameter.')
      endif
    case(ADJ_PARAM_ACTIVATION_ENERGY)
      microbial => realization%reaction%microbial
      if (microbial%nrxn == 0) then
        option%io_buffer = 'ACTIVATION_ENERGY KINETICS MICROBIAL specified, &
          &but the chemistry deck contains no MICROBIAL_REACTION block.'
        call PrintErrMsg(option)
      endif
      if (this%reaction_index < 1 .or. this%reaction_index > microbial%nrxn) &
        then
        option%io_buffer = 'REACTION_INDEX out of range for &
          &ACTIVATION_ENERGY KINETICS MICROBIAL.'
        call PrintErrMsg(option)
      endif
      ! the array exists when some reaction has one
      this%value = 0.d0
      if (associated(microbial%activation_energy)) &
        this%value = microbial%activation_energy(this%reaction_index)
      if (.not.ieee_is_finite(this%value) .or. this%value <= 0.d0) then
        call PrintErrMsg(option,'PARAMETER ACTIVATION_ENERGY requires a &
          &finite, positive ACTIVATION_ENERGY in the reaction''s &
          &MICROBIAL_REACTION block.')
      endif
    case(ADJ_PARAM_DECAY_RATE)
      if (realization%reaction%immobile%ndecay_rxn == 0) then
        option%io_buffer = 'DECAY_RATE_CONSTANT specified, but the &
          &chemistry deck contains no IMMOBILE_SPECIES / &
          &IMMOBILE_DECAY_REACTION block.'
        call PrintErrMsg(option)
      endif
      if (this%reaction_index < 1 .or. &
          this%reaction_index > realization%reaction%immobile%ndecay_rxn) then
        option%io_buffer = 'REACTION_INDEX out of range for &
          &DECAY_RATE_CONSTANT (must index an IMMOBILE_DECAY_REACTION).'
        call PrintErrMsg(option)
      endif
      this%value = &
        realization%reaction%immobile%decay_rate_constant(this%reaction_index)
      if (.not.ieee_is_finite(this%value) .or. this%value < 0.d0) then
        call PrintErrMsg(option,'DECAY adjoint requires a finite, &
          &nonnegative parameter.')
      endif
    case(ADJ_PARAM_PERMEABILITY)
      ! one isotropic value per material, at least one active cell
      mat => MaterialPropGetPtrFromList(this%material_name, &
                                        realization%material_properties)
      if (.not.associated(mat)) then
        call PrintErrMsg(option,'PARAMETER PERMEABILITY: material "' // &
                         trim(this%material_name) // '" not found.')
      endif
      if (associated(mat%permeability_dataset) .or. &
          associated(mat%permeability_dataset_y) .or. &
          associated(mat%permeability_dataset_z) .or. &
          associated(mat%permeability_dataset_xy) .or. &
          associated(mat%permeability_dataset_xz) .or. &
          associated(mat%permeability_dataset_yz)) then
        call PrintErrMsg(option,'PARAMETER PERMEABILITY does not support &
          &permeability datasets (one value per material).')
      endif
      if (.not.mat%isotropic_permeability) then
        call PrintErrMsg(option,'PARAMETER PERMEABILITY requires an &
          &isotropic material permeability (PERM_ISO).')
      endif
      this%material_id = mat%internal_id
      this%value = mat%permeability(1,1)
      if (.not.ieee_is_finite(this%value) .or. this%value <= 0.d0) then
        call PrintErrMsg(option,'PARAMETER PERMEABILITY requires a positive &
          &permeability.')
      endif
      patch => realization%patch
      this%value0 = this%value
      allocate(this%perm0(size(patch%aux%Material%auxvars(1)%permeability), &
                          patch%grid%ngmax))
      do local_id = 1, patch%grid%ngmax
        this%perm0(:,local_id) = &
          patch%aux%Material%auxvars(local_id)%permeability
      enddo
      ncells = 0
      do local_id = 1, patch%grid%nlmax
        if (patch%imat(patch%grid%nL2G(local_id)) == this%material_id) &
          ncells = ncells + 1
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,ncells,ONE_INTEGER_MPI,MPIU_INTEGER, &
                         MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
      if (ncells == 0) then
        call PrintErrMsg(option,'PARAMETER PERMEABILITY: no active cell is &
          &of material "' // trim(this%material_name) // '".')
      endif
    case(ADJ_PARAM_INITIAL_CONCENTRATION)
      reaction => realization%reaction
      id = ReactionAuxGetPriSpecIDFromName(this%species_name,reaction, &
                                           PETSC_FALSE,option)
      if (id > 0) then
        this%component = id
      else
        id = ReactionImGetSpeciesIDFromName(this%species_name, &
                                            reaction%immobile,PETSC_FALSE, &
                                            option)
        if (id <= 0) then
          call PrintErrMsg(option,'PARAMETER INITIAL_CONCENTRATION: "' // &
            trim(this%species_name) // '" is not a primary or immobile &
            &species.')
        endif
        this%component = reaction%offset_immobile + id
      endif
      patch => realization%patch
      region => RegionGetPtrFromList(this%region_name,patch%region_list)
      if (.not.associated(region)) then
        call PrintErrMsg(option,'PARAMETER INITIAL_CONCENTRATION: region "' &
          // trim(this%region_name) // '" not found.')
      endif
      ! the region's active cells here, and their initial value in the
      ! transport solution (free molality, immobile concentration)
      allocate(cells(region%num_cells))
      ncells = 0
      extremes = [huge(1.d0),huge(1.d0)]
      call VecGetArrayRead(realization%field%tran_xx,xx_p, &
                           ierr);CHKERRQ(ierr)
      do i = 1, region%num_cells
        local_id = region%cell_ids(i)
        if (patch%imat(patch%grid%nL2G(local_id)) <= 0) cycle
        ncells = ncells + 1
        cells(ncells) = local_id
        this%value = xx_p((local_id-1)*reaction%ncomp+this%component)
        ! min and -max
        extremes = min(extremes,[this%value,-this%value])
      enddo
      call VecRestoreArrayRead(realization%field%tran_xx,xx_p, &
                               ierr);CHKERRQ(ierr)
      this%cells = cells(1:ncells)
      call MPI_Allreduce(MPI_IN_PLACE,ncells,ONE_INTEGER_MPI,MPIU_INTEGER, &
                         MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
      call MPI_Allreduce(MPI_IN_PLACE,extremes,TWO_INTEGER_MPI, &
                         MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm, &
                         ierr);CHKERRQ(ierr)
      if (ncells == 0) then
        call PrintErrMsg(option,'PARAMETER INITIAL_CONCENTRATION: region "' &
          // trim(this%region_name) // '" has no active cell.')
      endif
      this%value = -extremes(2)
      if (.not.(extremes(1) > 0.d0 .and. ieee_is_finite(this%value))) then
        call PrintErrMsg(option,'PARAMETER ' // trim(this%name) // &
          ' requires a positive, finite initial concentration.')
      endif
      if (this%value > extremes(1)*(1.d0+1.d-12)) then
        write(option%io_buffer,'(a,es12.5,a,es12.5,a)') 'PARAMETER ' // &
          trim(this%name) // ': the initial concentrations in the region &
          &range from ',extremes(1),' to ',this%value,'; the parameter is &
          &one value.'
        call PrintErrMsg(option)
      endif
  end select

end subroutine AdjParamSetup

! ************************************************************************** !

subroutine AdjParamSelectTerm(this,realization,option)
  !
  ! The Monod (HALF_SATURATION_CONSTANT) or inhibition term of the
  ! parameter's reaction on its species: an error unless exactly one term
  ! matches. SMOOTHSTEP_INTERVAL considers SMOOTHSTEP terms only.
  !
  ! Author: Hansol Suh
  ! Date: 06/23/26
  !
  use Option_module
  use Realization_Subsurface_class
  use String_module, only : StringCompare
  use Reaction_Microbial_Aux_module, only : microbial_type
  use Reaction_Inhibition_Aux_module, only : INHIBITION_SMOOTHSTEP

  implicit none

  class(adjoint_param_type) :: this
  class(realization_subsurface_type) :: realization
  type(option_type) :: option

  type(microbial_type), pointer :: microbial
  PetscInt :: i, candidate, component, count, nterms, irxn
  PetscBool :: monod

  microbial => realization%reaction%microbial
  irxn = this%reaction_index
  monod = this%kind == ADJ_PARAM_HALF_SATURATION
  count = 0
  this%term = UNINITIALIZED_INTEGER
  if (monod) then
    nterms = microbial%monodid(0,irxn)
  else
    nterms = microbial%inhibitionid(0,irxn)
  endif
  do i = 1, nterms
    if (monod) then
      candidate = microbial%monodid(i,irxn)
      component = microbial%monod_specid(candidate)
    else
      candidate = microbial%inhibitionid(i,irxn)
      component = microbial%inhibition_specid(candidate)
      if (this%kind == ADJ_PARAM_SMOOTHSTEP_INTERVAL) then
        if (microbial%inhibition_type(candidate) /= INHIBITION_SMOOTHSTEP) &
          cycle
      endif
    endif
    if (.not.StringCompare(this%species_name, &
                           realization%reaction% &
                             primary_species_names(component), &
                           MAXWORDLENGTH)) cycle
    this%term = candidate
    count = count + 1
  enddo
  if (count > 1) then
    write(option%io_buffer,'(a,i0,a,a,a)') &
      'Ambiguous adjoint term: REACTION_INDEX ',irxn,', SPECIES_NAME "', &
      trim(this%species_name),'" matches multiple eligible blocks.'
    call PrintErrMsg(option)
  endif
  if (count == 0) then
    select case(this%kind)
      case(ADJ_PARAM_HALF_SATURATION)
        call PrintErrMsg(option,'HALF_SATURATION_CONSTANT: SPECIES_NAME &
          &does not match any MONOD block in the selected reaction.')
      case(ADJ_PARAM_SMOOTHSTEP_INTERVAL)
        call PrintErrMsg(option,'SMOOTHSTEP_INTERVAL: SPECIES_NAME does not &
          &match any TYPE SMOOTHSTEP inhibition block in the selected &
          &reaction.')
      case default
        call PrintErrMsg(option,'INHIBITION_CONSTANT: SPECIES_NAME does not &
          &match any INHIBITION block in the selected reaction.')
    end select
  endif

end subroutine AdjParamSelectTerm

! ************************************************************************** !

subroutine AdjParamSetValue(this,realization,val)
  !
  ! Sets value and the parameter in the reaction or material data the
  ! model reads; the caller invalidates the model's cache. A permeability
  ! scales the material's cell permeabilities, ghosts included, from their
  ! Setup values; the model rereads them in InvalidateCache. An initial
  ! concentration is only the value: the adjoint sets the initial state.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  use Realization_Subsurface_class
  use Patch_module

  implicit none

  class(adjoint_param_type) :: this
  class(realization_subsurface_type) :: realization
  PetscReal :: val

  type(patch_type), pointer :: patch
  PetscInt :: ghosted_id

  select case(this%kind)
    case(ADJ_PARAM_RATE_CONSTANT)
      if (this%kinetics == ADJ_KINETICS_MINERAL) then
        realization%reaction%mineral% &
          kinmnrl_precip_rate_constant(this%imnrl) = val
        realization%reaction%mineral% &
          kinmnrl_dissol_rate_constant(this%imnrl) = val
      else
        realization%reaction%microbial%rate_constant(this%reaction_index) = &
          val
      endif
    case(ADJ_PARAM_HALF_SATURATION)
      realization%reaction%microbial%monod_K(this%term) = val
    case(ADJ_PARAM_INHIBITION_CONSTANT)
      realization%reaction%microbial%inhibition_C(this%term) = &
        this%threshold_sign*val
    case(ADJ_PARAM_SMOOTHSTEP_INTERVAL)
      realization%reaction%microbial%inhibition_C2(this%term) = val
    case(ADJ_PARAM_YIELD)
      realization%reaction%microbial%biomass_yield(this%reaction_index) = val
    case(ADJ_PARAM_ACTIVATION_ENERGY)
      realization%reaction%microbial% &
        activation_energy(this%reaction_index) = val
    case(ADJ_PARAM_DECAY_RATE)
      realization%reaction%immobile% &
        decay_rate_constant(this%reaction_index) = val
    case(ADJ_PARAM_PERMEABILITY)
      patch => realization%patch
      do ghosted_id = 1, patch%grid%ngmax
        if (patch%imat(ghosted_id) /= this%material_id) cycle
        patch%aux%Material%auxvars(ghosted_id)%permeability = &
          (val/this%value0)*this%perm0(:,ghosted_id)
      enddo
  end select
  this%value = val

end subroutine AdjParamSetValue

! ************************************************************************** !

function AdjParamCount()
  !
  ! Returns the number of parameters currently in the list.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  implicit none

  PetscInt :: AdjParamCount
  type(adjoint_param_type), pointer :: cur

  AdjParamCount = 0
  cur => adj_param_list
  do
    if (.not.associated(cur)) exit
    AdjParamCount = AdjParamCount + 1
    cur => cur%next
  enddo

end function AdjParamCount

! ************************************************************************** !

subroutine AdjParamDestroyList()
  !
  ! Destroys every parameter in the list.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  implicit none

  type(adjoint_param_type), pointer :: cur, prev

  cur => adj_param_list
  do
    if (.not.associated(cur)) exit
    prev => cur%next
    deallocate(cur)
    cur => prev
  enddo
  nullify(adj_param_list)

end subroutine AdjParamDestroyList

end module Adjoint_Param_module
