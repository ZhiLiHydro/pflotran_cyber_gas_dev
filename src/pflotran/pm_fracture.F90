module PM_Fracture_class

#include "petsc/finclude/petscvec.h"
  use petscvec

  use Realization_Subsurface_class
  use PFLOTRAN_Constants_module
  use PM_Base_class
  use Option_module
  use String_module
  use Input_Aux_module
  use Geometry_module

  implicit none

  private

  type :: fracture_type
    ! fracture ID number
    PetscInt :: id
    ! cell permeability, initial [m]
    PetscReal, pointer :: kx0(:)
    ! cell permeability, initial [m]
    PetscReal, pointer :: ky0(:)
    ! cell permeability, initial [m]
    PetscReal, pointer :: kz0(:)
    ! fracture hydraulic aperture, initial [m]
    PetscReal :: hap0 ! 
    ! ejb : add dip here -- do Ctrl+F for hap0, can follow procedures for that
    type(point3d_type) :: strike ! direction along which the fracture plane intersects with the horizontal (assumed xz) plane
    PetscReal :: dip  ! angle in degrees between the horizontal plane and the fracture. goes along with strike
    ! fracture hydraulic aperture, dynamic [m]
    PetscReal, pointer :: hap(:)
    ! fracture radius - x [m]
    PetscReal :: radx
    ! fracture radius - y [m]
    PetscReal :: rady
    ! fracture radius - z [m]
    PetscReal :: radz
    ! fracture center coordinates [m]
    type(point3d_type) :: center
    ! fracture normal vector [m]
    type(point3d_type) :: normal
    ! fracture parallel vector [m]
    type(point3d_type) :: parallel
    ! rotation matrix
    PetscReal :: RM(3,3)
    ! list of local cell ids the fracture occupies
    PetscInt, pointer :: cell_ids(:)
    ! list of previous temperatures in the fracture cells
    PetscReal, pointer :: prev_temperature(:)
    ! list of change in temperatures in the fracture cells
    PetscReal, pointer :: dT(:)
    ! list of change in length in the fracture cells
    PetscReal, pointer :: dL(:)
    ! total number of cells this fracture occupies
    PetscInt :: ncells
    ! linked list next object
    type(fracture_type), pointer :: next ! ejb: pointer to the next fracture -- know you hit the end when the next pointer is unassociated
  end type fracture_type

  type :: fracfam_type
    ! fracture family ID number
    PetscInt :: id
    ! number of fractures in the family
    PetscInt :: nfrac_in_fam
    ! fracture intensity of fracture family [num fractures]
    PetscInt :: intensity_fracfam
    ! total surface area of fracture family [m^2] (includes both sides of fracture)
    PetscReal :: surface_area_fracfam
    ! fracture family hydraulic aperture, initial [m]
    PetscReal :: hap0
    ! fracture family hydraulic aperture, initial maximum [m]
    PetscReal :: hap0_max
    ! fracture family hydraulic standard deviation, initial [m]
    PetscReal :: hap0_stdev
    ! fracture family hydraulic aperture seed [-]
    PetscInt :: hap0_seed
    ! fracture family center coordinates [m]
    type(point3d_type) :: center
    ! fracture family center coordinates standard deviation [m]
    type(point3d_type) :: center_stdev
    ! fracture family center seed [-]
    PetscInt :: center_seed
    ! fracture family normal vector [m]
    type(point3d_type) :: normal
    ! fracture family normal vector coordinates standard deviation [m]
    type(point3d_type) :: normal_stdev
    ! fracture family normal vector seed [-]
    PetscInt :: normal_seed
    ! fracture family radius lengths [m]
    type(point3d_type) :: radius
    ! fracture family radius lengths standard deviation [m]
    type(point3d_type) :: radius_stdev
    ! fracture family radius seed [-]
    PetscInt :: radius_seed
    ! list of fractures in this family
    type(fracture_type), pointer :: fracture_list 
    ! linked list next object
    type(fracfam_type), pointer :: next
  end type fracfam_type

  type, public, extends(pm_base_type) :: pm_fracture_type
    class(realization_subsurface_type), pointer :: realization
    type(fracture_type), pointer :: fracture_list ! ejb:list of all the fractures in the domain (this includes fractures created within a fracture family) (call it global fracture list)
    type(fracfam_type), pointer :: fracfam_list
    ! list of all local cell ids that contain fractures (some repeated)
    PetscInt, pointer :: allfrac_cell_ids(:)
    ! list of all local cell ids that contain fractures (unique)
    PetscInt, pointer :: allfrac_unique_cell_ids(:)
    ! list of all global cell ids that contain fractures (unique)
    PetscInt, pointer :: allfrac_unique_cell_ids_glb(:)
    ! list of all local cell ids that contain fracture intersections
    PetscInt, pointer :: frac_common_cell_ids(:)
    ! Flag indicating connection to starting region
    PetscBool, pointer :: connected_to_start(:)
    ! Flag indicating connection to ending region
    PetscBool, pointer :: connected_to_end(:)
    ! list of all local cell ids that contain connected fractures (unique)
    PetscInt, pointer :: unique_connected_cell_ids(:)
    ! list of all global cell ids that contain connected fractures (unique)
    PetscInt, pointer :: unique_connected_cell_ids_glb(:)
    ! list of start region global cell ids
    PetscInt, pointer :: clip_region_cell_ids(:)
    ! name of region for clipping
    character(len=MAXWORDLENGTH) :: start_region_name
    ! name of end region for clipping
    character(len=MAXWORDLENGTH) :: end_region_name
    ! maximum fracture perm at global cell ids
    PetscReal, pointer :: max_frac_kx(:)
    PetscReal, pointer :: max_frac_ky(:)
    PetscReal, pointer :: max_frac_kz(:)
    ! sum of fracture perm at global cell ids
    PetscReal, pointer :: sum_frac_kx(:)
    PetscReal, pointer :: sum_frac_ky(:)
    PetscReal, pointer :: sum_frac_kz(:)
    ! maximum fracture aperture at global cell ids
    PetscReal, pointer :: max_frac_hap(:)
    ! sum of fracture aperture at global cell ids
    PetscReal, pointer :: sum_frac_hap(:)
    ! global cell ids where fracture intersections occur
    PetscInt, pointer :: frac_intersection_ids(:)
    ! behavior at fracture intersections (1=sum,2=max)
    PetscInt :: frac_intersection_type
    ! thermal expansion coefficient # [1/C]
    PetscReal :: t_coeff
    ! number of fractures in the domain
    PetscInt :: nfrac
    ! maximum number of fractures in the domain
    PetscInt :: max_frac
    ! flag to update material_auxvar (permeability) object
    PetscBool :: update_material_auxvar_perm
    PetscBool :: clip_unconnected_fracs
  contains
    procedure, public :: Setup => PMFracSetup
    procedure, public :: ReadPMBlock => PMFracReadPMBlock
    procedure, public :: SetRealization => PMFracSetRealization
    procedure, public :: InitializeRun => PMFracInitializeRun
    procedure, public :: FinalizeRun => PMFracFinalizeRun
    procedure, public :: InitializeTimestep => PMFracInitializeTimestep
    !procedure, public :: UpdateTimestep => PMFracUpdateTimestep
    procedure, public :: FinalizeTimestep => PMFracFinalizeTimestep
    !procedure, public :: PreSolve => PMFracPreSolve
    procedure, public :: Solve => PMFracSolve
    !procedure, public :: PostSolve => PMFracPostSolve
    procedure, public :: Destroy => PMFracDestroy
  end type pm_fracture_type

  public :: PMFractureCreate, &
            FractureCreate, &
            FractureFamCreate, &
            PMFracReadPMBlock, &
            PMFracReadPass2

  contains

! ************************************************************************** !

function PMFractureCreate()
  !
  ! Creates the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023

  implicit none

  class(pm_fracture_type), pointer :: PMFractureCreate
  class(pm_fracture_type), pointer :: this

  allocate(this)
  call PMBaseInit(this)

  this%header = 'GEOTHERMAL FRACTURE MODEL'

  nullify(this%realization)
  nullify(this%fracture_list)
  nullify(this%fracfam_list)
  nullify(this%allfrac_cell_ids)
  nullify(this%allfrac_unique_cell_ids)
  nullify(this%allfrac_unique_cell_ids_glb)
  nullify(this%frac_common_cell_ids)
  nullify(this%connected_to_start)
  nullify(this%connected_to_end)
  nullify(this%unique_connected_cell_ids)
  nullify(this%unique_connected_cell_ids_glb)
  nullify(this%clip_region_cell_ids)
  nullify(this%max_frac_kx)
  nullify(this%max_frac_ky)
  nullify(this%max_frac_kz)
  nullify(this%sum_frac_kx)
  nullify(this%sum_frac_ky)
  nullify(this%sum_frac_kz)
  nullify(this%max_frac_hap)
  nullify(this%sum_frac_hap)
  nullify(this%frac_intersection_ids)
  this%frac_intersection_type = 1  ! (1=sum;2=max)
  this%nfrac = 0
  this%max_frac = UNINITIALIZED_INTEGER
  this%t_coeff = UNINITIALIZED_DOUBLE ! 40.d-6 is value for granite
  this%update_material_auxvar_perm = PETSC_FALSE
  this%clip_unconnected_fracs = PETSC_FALSE
  this%start_region_name = ''
  this%end_region_name = ''

  PMFractureCreate => this

end function PMFractureCreate

! ************************************************************************** !

function FractureCreate()
  !
  ! Creates a fracture object.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023

  implicit none

  type(fracture_type), pointer :: FractureCreate

  allocate(FractureCreate)
  call FractureInit(FractureCreate)

end function FractureCreate

! ************************************************************************** !

function FractureFamCreate()
  !
  ! Creates a fracture family object.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/19/2024

  implicit none

  type(fracfam_type), pointer :: FractureFamCreate

  allocate(FractureFamCreate)
  call FractureFamInit(FractureFamCreate)

end function FractureFamCreate

! ************************************************************************** !

subroutine FractureInit(this)
  !
  ! Initializes variables associated with a fracture object.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !
  implicit none

  type(fracture_type) :: this

  nullify(this%next)
  nullify(this%cell_ids)
  nullify(this%prev_temperature)
  nullify(this%dT)
  nullify(this%dL)
  nullify(this%hap)
  nullify(this%kx0); nullify(this%ky0); nullify(this%kz0)
  this%ncells = UNINITIALIZED_INTEGER
  this%id = UNINITIALIZED_INTEGER
  this%hap0 = UNINITIALIZED_DOUBLE
  this%dip = UNINITIALIZED_DOUBLE
  this%radx = UNINITIALIZED_DOUBLE
  this%rady = UNINITIALIZED_DOUBLE
  this%radz = UNINITIALIZED_DOUBLE
  this%center%x = UNINITIALIZED_DOUBLE
  this%center%y = UNINITIALIZED_DOUBLE
  this%center%z = UNINITIALIZED_DOUBLE
  this%normal%x = UNINITIALIZED_DOUBLE
  this%normal%y = UNINITIALIZED_DOUBLE
  this%normal%z = UNINITIALIZED_DOUBLE
  this%strike%x = UNINITIALIZED_DOUBLE
  this%strike%y = UNINITIALIZED_DOUBLE
  this%strike%z = UNINITIALIZED_DOUBLE
  this%parallel%x = UNINITIALIZED_DOUBLE
  this%parallel%y = UNINITIALIZED_DOUBLE
  this%parallel%z = UNINITIALIZED_DOUBLE

end subroutine FractureInit

! ************************************************************************** !

subroutine FractureFamInit(this)
  !
  ! Initializes variables associated with a fracture family object.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/19/2024
  !
  implicit none

  type(fracfam_type) :: this

  nullify(this%next)
  nullify(this%fracture_list)
  this%id = UNINITIALIZED_INTEGER
  this%nfrac_in_fam = 0
  this%intensity_fracfam= UNINITIALIZED_DOUBLE
  this%surface_area_fracfam = 0.d0
  this%hap0 = UNINITIALIZED_DOUBLE
  this%hap0_stdev = UNINITIALIZED_DOUBLE
  this%hap0_seed = 5
  this%hap0_max = UNINITIALIZED_DOUBLE
  this%center%x = UNINITIALIZED_DOUBLE
  this%center%y = UNINITIALIZED_DOUBLE
  this%center%z = UNINITIALIZED_DOUBLE
  this%center_stdev%x = UNINITIALIZED_DOUBLE
  this%center_stdev%y = UNINITIALIZED_DOUBLE
  this%center_stdev%z = UNINITIALIZED_DOUBLE
  this%center_seed = 10
  this%normal%x = UNINITIALIZED_DOUBLE
  this%normal%y = UNINITIALIZED_DOUBLE
  this%normal%z = UNINITIALIZED_DOUBLE
  this%normal_stdev%x = UNINITIALIZED_DOUBLE
  this%normal_stdev%y = UNINITIALIZED_DOUBLE
  this%normal_stdev%z = UNINITIALIZED_DOUBLE
  this%normal_seed = 20
  this%radius%x = UNINITIALIZED_DOUBLE
  this%radius%y = UNINITIALIZED_DOUBLE
  this%radius%z = UNINITIALIZED_DOUBLE
  this%radius_stdev%x = UNINITIALIZED_DOUBLE
  this%radius_stdev%y = UNINITIALIZED_DOUBLE
  this%radius_stdev%z = UNINITIALIZED_DOUBLE
  this%radius_seed = 30

end subroutine FractureFamInit

! ************************************************************************** !

subroutine PMFracSetup(this)
  !
  ! Initializes variables associated with the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !
  use Region_module
  use Material_Aux_module
  use Grid_module

  implicit none

  class(pm_fracture_type) :: this

  type(option_type), pointer :: option
  type(fracture_type), pointer :: cur_fracture,new_fracture_list
  type(fracfam_type), pointer :: cur_fracfam
  type(grid_type), pointer :: res_grid
  type(material_auxvar_type), pointer :: material_auxvar
  character(len=MAXWORDLENGTH) :: word,filename, filename_frac_info
  character(len=MAXSTRINGLENGTH) :: path
  PetscInt, pointer :: temp_cell_ids(:),temp_gbl_cell_ids(:)
  PetscInt, pointer :: temp_allfrac_cell_ids(:)
  PetscInt, pointer :: temp_allfrac_unique_cell_ids(:)
  PetscInt, pointer :: temp_frac_common_cell_ids(:)
  PetscInt :: k,j,i,fid,status, fid_frac_info, glb_id
  PetscInt :: nf,tfc,tcfc,read_max
  PetscInt :: total_sum_glb
  PetscInt, pointer :: allfrac_num_cells(:),common_num_cells(:)
  PetscInt, pointer :: uniq_conn_num_cells(:)
  PetscInt, pointer :: cell_positions(:)
  PetscReal, pointer :: d_vertex(:)
  PetscReal :: D,distance,L22
  PetscReal :: min_x,max_x,min_y,max_y,min_z,max_z
  PetscReal :: a1,a2,a3,b1,b2,b3,c1,c2,c3,vmag
  PetscReal :: ra, sinra, cosra
  type(point3d_type) :: grid_center,grid_dh ! grid cell center point and dh
  type(point3d_type) :: v1,v2,v3,v4,v5,v6,v7,v8 ! grid cell vertices
  PetscBool :: within_x,within_y,within_z
  PetscBool :: frac_is_within_grid_cell_volume
  PetscBool :: added
  PetscBool, pointer :: connected(:)
  PetscErrorCode :: ierr

  call this%SetRealization()

  option => this%option
  res_grid => this%realization%patch%grid

  do k = 1,res_grid%nlmax
    material_auxvar => &
        this%realization%patch%aux%material%auxvars(res_grid%nL2G(k))
    material_auxvar%fracture_aperture = 0.d0
  enddo

  nf = 0
  read_max = 15000 ! Set to a large but realistic number

  allocate(temp_cell_ids(read_max))
  allocate(temp_allfrac_cell_ids(read_max*10))
  allocate(temp_frac_common_cell_ids(read_max))
  allocate(d_vertex(8))

  cur_fracfam => this%fracfam_list
  do
    if (.not.associated(cur_fracfam)) exit

    ! generate fractures in this fracture family
    call PMFracGenerateFracFam(cur_fracfam)

    cur_fracfam => cur_fracfam%next
  enddo

  added = PETSC_FALSE
  cur_fracfam => this%fracfam_list
  do
    if (.not.associated(cur_fracfam)) exit
    new_fracture_list => cur_fracfam%fracture_list

    ! add fracture family fractures to the individual fracture list
    if (.not.associated(this%fracture_list)) then
      this%fracture_list => new_fracture_list
    else
      cur_fracture => this%fracture_list
      do
        if (.not.associated(cur_fracture)) exit
        if (.not.associated(cur_fracture%next)) then
          cur_fracture%next => new_fracture_list
          added = PETSC_TRUE
        endif
        if (added) exit
        cur_fracture => cur_fracture%next
      enddo
    endif

    cur_fracfam => cur_fracfam%next
    added = PETSC_FALSE
  enddo

  cur_fracture => this%fracture_list
  do
    if (.not.associated(cur_fracture)) exit
    nf = nf + 1
    cur_fracture => cur_fracture%next
  enddo
  write(word,'(i4)') nf

  option%io_buffer = ' '
  call PrintMsg(option)
  option%io_buffer = 'GEOTHERMAL_FRACTURE_MODEL: Mapping [' // trim(word) // &
    '] fractures onto equivalent porous media domain:'
  call PrintMsg(option)

  tfc = 0
  cur_fracture => this%fracture_list
  do
    if (.not.associated(cur_fracture)) exit

    write(word,'(i4)') cur_fracture%id
    option%io_buffer = 'GEOTHERMAL_FRACTURE_MODEL: Mapping fracture ID# [' &
                       // trim(word) // '].'
    call PrintMsg(option)

    ! Equation of a plane normal to vector (A,B,C) is
    ! Ax + By + Cz + D = 0
    ! Get D by knowing the plane must contain the center point (x,y,z)
    D = -1.d0*(cur_fracture%normal%x*cur_fracture%center%x + &
               cur_fracture%normal%y*cur_fracture%center%y + &
               cur_fracture%normal%z*cur_fracture%center%z)
    j = 0
    do k = 1,res_grid%nlmax
      grid_center%x = res_grid%x(res_grid%nL2G(k)) ! grid cell center point x
      grid_center%y = res_grid%y(res_grid%nL2G(k)) ! grid cell center point y
      grid_center%z = res_grid%z(res_grid%nL2G(k)) ! grid cell center point z
      if (res_grid%itype == STRUCTURED_GRID) then
        grid_dh%x = res_grid%structured_grid%dx(res_grid%nL2G(k))
        grid_dh%y = res_grid%structured_grid%dy(res_grid%nL2G(k))
        grid_dh%z = res_grid%structured_grid%dz(res_grid%nL2G(k))
        ! calculate cuboid vertex points
        v1%x = grid_center%x - (grid_dh%x/2.d0)
        v2%x = grid_center%x + (grid_dh%x/2.d0)
        v3%x = grid_center%x - (grid_dh%x/2.d0)
        v4%x = grid_center%x + (grid_dh%x/2.d0)
        v5%x = grid_center%x - (grid_dh%x/2.d0)
        v6%x = grid_center%x + (grid_dh%x/2.d0)
        v7%x = grid_center%x - (grid_dh%x/2.d0)
        v8%x = grid_center%x + (grid_dh%x/2.d0)
        v1%y = grid_center%y - (grid_dh%y/2.d0)
        v2%y = grid_center%y - (grid_dh%y/2.d0)
        v3%y = grid_center%y + (grid_dh%y/2.d0)
        v4%y = grid_center%y + (grid_dh%y/2.d0)
        v5%y = grid_center%y - (grid_dh%y/2.d0)
        v6%y = grid_center%y - (grid_dh%y/2.d0)
        v7%y = grid_center%y + (grid_dh%y/2.d0)
        v8%y = grid_center%y + (grid_dh%y/2.d0)
        v1%z = grid_center%z - (grid_dh%z/2.d0)
        v2%z = grid_center%z - (grid_dh%z/2.d0)
        v3%z = grid_center%z - (grid_dh%z/2.d0)
        v4%z = grid_center%z - (grid_dh%z/2.d0)
        v5%z = grid_center%z + (grid_dh%z/2.d0)
        v6%z = grid_center%z + (grid_dh%z/2.d0)
        v7%z = grid_center%z + (grid_dh%z/2.d0)
        v8%z = grid_center%z + (grid_dh%z/2.d0)
      else
        option%io_buffer = 'The GEOTHERMAL_FRACTURE_MODEL can only be used &
          &with structured grids. If you recieve this error message. . .'
        call PrintErrMsgToDev(option, &
          'describe your modeling application or goals.')
      endif
      d_vertex(1) = cur_fracture%normal%x*v1%x + &
                    cur_fracture%normal%y*v1%y + &
                    cur_fracture%normal%z*v1%z + D
      d_vertex(2) = cur_fracture%normal%x*v2%x + &
                    cur_fracture%normal%y*v2%y + &
                    cur_fracture%normal%z*v2%z + D
      d_vertex(3) = cur_fracture%normal%x*v3%x + &
                    cur_fracture%normal%y*v3%y + &
                    cur_fracture%normal%z*v3%z + D
      d_vertex(4) = cur_fracture%normal%x*v4%x + &
                    cur_fracture%normal%y*v4%y + &
                    cur_fracture%normal%z*v4%z + D
      d_vertex(5) = cur_fracture%normal%x*v5%x + &
                    cur_fracture%normal%y*v5%y + &
                    cur_fracture%normal%z*v5%z + D
      d_vertex(6) = cur_fracture%normal%x*v6%x + &
                    cur_fracture%normal%y*v6%y + &
                    cur_fracture%normal%z*v6%z + D
      d_vertex(7) = cur_fracture%normal%x*v7%x + &
                    cur_fracture%normal%y*v7%y + &
                    cur_fracture%normal%z*v7%z + D
      d_vertex(8) = cur_fracture%normal%x*v8%x + &
                    cur_fracture%normal%y*v8%y + &
                    cur_fracture%normal%z*v8%z + D
      distance = cur_fracture%normal%x*grid_center%x + &
                 cur_fracture%normal%y*grid_center%y + &
                 cur_fracture%normal%z*grid_center%z + D
      frac_is_within_grid_cell_volume = PETSC_TRUE
      if ((all(d_vertex > 0.d0)) .or. (all(d_vertex < 0.d0))) then
        ! if the distances between all vertices and the fracture plane are
        ! all positive or all negative, then the grid cell volume is on either
        ! side of the plane, and the fracture plane doesn't go through
        frac_is_within_grid_cell_volume = PETSC_FALSE
      elseif (any(d_vertex == 0.d0)) then
        ! if any vertex is exactly on the plane. . . .
        frac_is_within_grid_cell_volume = PETSC_TRUE
      endif
      if (frac_is_within_grid_cell_volume) then
        min_x = cur_fracture%center%x - cur_fracture%radx
        max_x = cur_fracture%center%x + cur_fracture%radx
        min_y = cur_fracture%center%y - cur_fracture%rady
        max_y = cur_fracture%center%y + cur_fracture%rady
        min_z = cur_fracture%center%z - cur_fracture%radz
        max_z = cur_fracture%center%z + cur_fracture%radz
        within_x = PETSC_FALSE; within_y = PETSC_FALSE
        within_z = PETSC_FALSE
        if ((grid_center%x <= max_x) .and. (grid_center%x >= min_x)) then
          within_x = PETSC_TRUE
        endif
        if ((grid_center%y <= max_y) .and. (grid_center%y >= min_y)) then
          within_y = PETSC_TRUE
        endif
        if ((grid_center%z <= max_z) .and. (grid_center%z >= min_z)) then
          within_z = PETSC_TRUE
        endif
        ! if you want the fracture oriented exactly along an axis, then allow
        ! that axis to pass the following logic statement by setting the
        ! within_* booleans to TRUE
        if (abs(max_x-min_x) < 1.d-5) within_x = PETSC_TRUE
        if (abs(max_y-min_y) < 1.d-5) within_y = PETSC_TRUE
        if (abs(max_z-min_z) < 1.d-5) within_z = PETSC_TRUE
        if (within_x .and. within_y .and. within_z) then
          j = j + 1
          temp_cell_ids(j) = k ! local cell id
        endif
      endif
    enddo
    cur_fracture%ncells = j
    if (j > read_max) then
      option%io_buffer = 'The number of grid cells that FRACTURE with &
        &ID# [' // word // '] occupies exceeds the set limit of 15,000 cells. &
        &If you are certain the fracture should be this large, increase the &
        &size of read_max in PMFracSetup() or email the PFLOTRAN developers &
        &group.'
      call PrintErrMsg(option)
    endif
    if (cur_fracture%ncells > 0) then
      allocate(cur_fracture%cell_ids(cur_fracture%ncells))
      allocate(cur_fracture%prev_temperature(cur_fracture%ncells))
      allocate(cur_fracture%dT(cur_fracture%ncells))
      allocate(cur_fracture%dL(cur_fracture%ncells))
      allocate(cur_fracture%hap(cur_fracture%ncells))
      allocate(cur_fracture%kx0(cur_fracture%ncells))
      allocate(cur_fracture%ky0(cur_fracture%ncells))
      allocate(cur_fracture%kz0(cur_fracture%ncells))

      cur_fracture%hap(:) = cur_fracture%hap0
      cur_fracture%cell_ids(1:cur_fracture%ncells) = &
                                          temp_cell_ids(1:cur_fracture%ncells)
      temp_allfrac_cell_ids(tfc+1:tfc+cur_fracture%ncells) = &
                                  cur_fracture%cell_ids(1:cur_fracture%ncells)
    endif

    tfc = tfc + cur_fracture%ncells
    if (tfc > (read_max*10)) then
      option%io_buffer = 'tfc exceed maximum. Email the PFLOTRAN developers &
        &group, or increase the size of read_max*10 in PMFracSetup().'
      call PrintErrMsg(option)
    endif

    ! get the unit vector parallel to the fracture plane
    ! a X b = c, where c is perpendicular to a and b (a and b can't be parallel)
    ! take a as the unit normal, and b as slightly rotated unit normal
    ! then c will be perpendicular to a, and parallel to fracture plane
    vmag = sqrt((cur_fracture%normal%x)**2.d0 + &
                (cur_fracture%normal%y)**2.d0 + &
                (cur_fracture%normal%z)**2.d0)
    cur_fracture%normal%x = cur_fracture%normal%x/vmag
    cur_fracture%normal%y = cur_fracture%normal%y/vmag
    cur_fracture%normal%z = cur_fracture%normal%z/vmag
    a1 = cur_fracture%normal%x; b1 = 0.d0
    a2 = cur_fracture%normal%y; b2 = 0.d0
    a3 = cur_fracture%normal%z; b3 = 1.d0
    ! take the cross product, c = a x b
    c1 = a2*b3 - a3*b2
    c2 = a3*b1 - a1*b3
    c3 = a1*b2 - a2*b1

    cur_fracture%parallel%x = c1
    cur_fracture%parallel%y = c2
    cur_fracture%parallel%z = c3
    vmag = sqrt((cur_fracture%parallel%x)**2.d0 + &
              (cur_fracture%parallel%y)**2.d0 + &
              (cur_fracture%parallel%z)**2.d0)
    if (vmag > 1.d-20) then
      cur_fracture%parallel%x = cur_fracture%parallel%x/vmag
      cur_fracture%parallel%y = cur_fracture%parallel%y/vmag
      cur_fracture%parallel%z = cur_fracture%parallel%z/vmag
    endif
    c1 = cur_fracture%parallel%x
    c2 = cur_fracture%parallel%y
    c3 = cur_fracture%parallel%z

    ! calculate the rotation angle (ra) between domain/fracture orientations
    ra = acos(a3)

    ! calculate the rotation matrix between domain/fracture orientations
    ! this is done knowing the rotation angle, and fracture parallel direction
    sinra = sin(ra); cosra = cos(ra)
    cur_fracture%RM(1,1) = cosra + c1*c1*(1.d0-cosra)
    cur_fracture%RM(1,2) = c1*c2*(1.d0-cosra) - c3*sinra
    cur_fracture%RM(1,3) = c1*c3*(1.d0-cosra) + c2*sinra
    cur_fracture%RM(2,1) = c2*c1*(1.d0-cosra) + c3*sinra
    cur_fracture%RM(2,2) = cosra + c2*c2*(1.d0-cosra)
    cur_fracture%RM(2,3) = c2*c3*(1.d0-cosra) -c1*sinra
    cur_fracture%RM(3,1) = c3*c1*(1.d0-cosra) -c2*sinra
    cur_fracture%RM(3,2) = c3*c2*(1.d0-cosra) + c1*sinra
    cur_fracture%RM(3,3) = cosra + c3*c3*(1.d0-cosra)

    ! calculate the dip and strike of the fracture 
    ! where strike is: the direction along which the fracture plane intersects the horizontal plane
    ! and dip is: the angle the fracture plane makes with the horizontal plane
    ! in our case, we assume the horizontal plane is the xy plane, with unit normal n = <0, 0, 1>
    ! Letting nfrac = unit vector for the fracture, nxy = xy plane unit vector, and nstrike = unit vector in strike direction, we have: 
    ! nstrike = nfrac X nxy =<(nfrac_y*1), -(nfrac_x*1), (0)> = <nfrac_y, -nfrac_x>
    cur_fracture%strike%x = cur_fracture%normal%y
    cur_fracture%strike%y = -1.d0*cur_fracture%normal%y
    cur_fracture%strike%z = 0.d0

    ! Then for dip, we have : ||nfrac||||nxy||cos(theta) = nfrac dot nxy, where theta is the dip angle
    ! since ||nfrac|| = ||nxy|| = 1, and nxy = <0,0,1> we have
    cur_fracture%dip = acos(cur_fracture%normal%z)


    cur_fracture => cur_fracture%next
  enddo

  allocate(this%allfrac_cell_ids(tfc)) ! list of local cell ids, ok if tfc=0
  this%allfrac_cell_ids(:) = UNINITIALIZED_INTEGER
  if (tfc > 0) this%allfrac_cell_ids(1:tfc) = temp_allfrac_cell_ids(1:tfc)
  deallocate(temp_allfrac_cell_ids)
  allocate(temp_allfrac_cell_ids(tfc)) ! reuse the array, now re-sized
  temp_allfrac_cell_ids(:) = UNINITIALIZED_INTEGER
  deallocate(temp_cell_ids)
  allocate(temp_cell_ids(tfc)) ! reuse the array, now re-sized
  temp_cell_ids(:) = UNINITIALIZED_INTEGER
  if (tfc > 0) then
    ! sort list of local cell ids in this%allfrac_cell_ids
    do i = 2,tfc
      k = this%allfrac_cell_ids(i)
      j = i - 1
      do while (j >= 1 .and. this%allfrac_cell_ids(j) > k)
        this%allfrac_cell_ids(j+1) = this%allfrac_cell_ids(j)
        j = j - 1
      end do
      this%allfrac_cell_ids(j+1) = k
    end do
    ! pull out the duplicated/repeated cell ids (these are intersections)
    k = 0; i = 1
    do while (i <= tfc)
      j = 1
      do while (((i + j) <= tfc) .and. &
                (this%allfrac_cell_ids(i+j) == this%allfrac_cell_ids(i)))
        j = j + 1
      end do
      if (j >= 2) then
        k = k + 1
        temp_cell_ids(k) = this%allfrac_cell_ids(i)
      end if
      i = i + j
    end do
    allocate(this%frac_common_cell_ids(k))
    if (k > 0) this%frac_common_cell_ids = temp_cell_ids(1:k)
  endif
  deallocate(temp_frac_common_cell_ids)

  ! if in parallel, the list in this%frac_common_cell_ids are local to the
  ! processor, so we need to make a global id list from all processors
  ! j is the total number of common cells on each process (ie intersections)
  allocate(common_num_cells(option%comm%size))
  j = size(this%frac_common_cell_ids)
  allocate(temp_gbl_cell_ids(j)) ! ok if j=0
  do i=1,j ! will not enter if j=0
    temp_gbl_cell_ids(i) = &
                     res_grid%nG2A(res_grid%nL2G(this%frac_common_cell_ids(i)))
  enddo
  call MPI_Allgather(j,ONE_INTEGER,MPI_INTEGER,common_num_cells, &
                     ONE_INTEGER,MPI_INTEGER,option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Barrier(option%mycomm,ierr);CHKERRQ(ierr)
  allocate(cell_positions(option%comm%size))
  do i = 1,option%comm%size
    if (i == 1) then
      cell_positions(i) = 0
    else
      cell_positions(i) = cell_positions(i-1) + common_num_cells(i-1)
    endif
  enddo
  allocate(this%frac_intersection_ids(SUM(common_num_cells))) ! ok if SUM=0
  if (SUM(common_num_cells) > 0) then
    call MPI_Gatherv(temp_gbl_cell_ids,j,MPI_INTEGER, &
                     this%frac_intersection_ids,common_num_cells, &
                     cell_positions,MPI_INTEGER,ZERO_INTEGER, &
                     option%mycomm,ierr);CHKERRQ(ierr)
    ! only rank zero knows this%frac_intersection_ids, so we need to share it
    call MPI_Barrier(option%mycomm,ierr);CHKERRQ(ierr)
    call MPI_Bcast(this%frac_intersection_ids,SUM(common_num_cells), &
                   MPI_INTEGER,ZERO_INTEGER,option%mycomm,ierr);CHKERRQ(ierr)
  endif
  deallocate(common_num_cells)
  deallocate(temp_gbl_cell_ids)

  allocate(temp_allfrac_unique_cell_ids(tfc)) ! list of local cell ids
  temp_allfrac_unique_cell_ids(:) = UNINITIALIZED_INTEGER
  j = 0
  if (tfc > 0) then
    temp_allfrac_unique_cell_ids(1) = this%allfrac_cell_ids(1)
    j = j + 1
    do i=2,tfc
      ! if the number already exists in temp check next
      if (any( temp_allfrac_unique_cell_ids == &
               this%allfrac_cell_ids(i) )) cycle
      ! No match found so add it to the output
      j = j + 1
      temp_allfrac_unique_cell_ids(j) = this%allfrac_cell_ids(i)
    enddo
  endif
  allocate(this%allfrac_unique_cell_ids(j)) ! ok if j=0
  this%allfrac_unique_cell_ids(:) = UNINITIALIZED_INTEGER
  if (j > 0) then
    this%allfrac_unique_cell_ids(1:j) = temp_allfrac_unique_cell_ids(1:j)
  endif
  ! reuse already allocated array to hold global ids before deallocating it
  do i=1,j ! will not enter if j=0
    temp_allfrac_unique_cell_ids(i) = &
                  res_grid%nG2A(res_grid%nL2G(this%allfrac_unique_cell_ids(i)))
  enddo
  ! if in parallel, the list in this%allfrac_unique_cell_ids are local to the
  ! processor, so we need to make a global id list from all processors
  ! j is the total number of fractured cells on each process
  allocate(allfrac_num_cells(option%comm%size))
  call MPI_Allgather(j,ONE_INTEGER,MPI_INTEGER,allfrac_num_cells, &
                     ONE_INTEGER,MPI_INTEGER,option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Barrier(option%mycomm,ierr);CHKERRQ(ierr)
  total_sum_glb = SUM(allfrac_num_cells)
  do i = 1,option%comm%size
    if (i == 1) then
      cell_positions(i) = 0
    else
      cell_positions(i) = cell_positions(i-1) + allfrac_num_cells(i-1)
    endif
  enddo
  allocate(this%allfrac_unique_cell_ids_glb(total_sum_glb)) ! ok if SUM=0
  if (total_sum_glb > 0) then
    call MPI_Gatherv(temp_allfrac_unique_cell_ids,j,MPI_INTEGER, &
                     this%allfrac_unique_cell_ids_glb,allfrac_num_cells, &
                     cell_positions,MPI_INTEGER,ZERO_INTEGER, &
                     option%mycomm,ierr);CHKERRQ(ierr)
    ! only rank zero knows this%allfrac_unique_cell_ids_glb, so we need to share it
    call MPI_Barrier(option%mycomm,ierr);CHKERRQ(ierr)
    call MPI_Bcast(this%allfrac_unique_cell_ids_glb,total_sum_glb,MPI_INTEGER, &
                   ZERO_INTEGER,option%mycomm,ierr);CHKERRQ(ierr)
  endif
  deallocate(allfrac_num_cells)
  deallocate(temp_allfrac_unique_cell_ids)

  if (this%clip_unconnected_fracs) then
    option%io_buffer = 'GEOTHERMAL_FRACTURE_MODEL: Calculating connected &
                       &fractures . . . '
    call PrintMsg(option)
    if (option%comm%size > 1) then
      option%io_buffer = 'GEOTHERMAL_FRACTURE_MODEL: The feature to clip &
        &unconnected fractures is not yet available when running in parallel.'
      call PrintErrMsg(option)
    endif

    !start region
    call PMFracClipCells(this,this%connected_to_start,total_sum_glb, &
                         cell_positions,this%start_region_name)
    !end region
    call PMFracClipCells(this,this%connected_to_end,total_sum_glb, &
                         cell_positions,this%end_region_name)

    allocate(connected(total_sum_glb))

    connected = this%connected_to_start .and. this%connected_to_end








    this%connected_to_start = connected









    !write(*,*) "............................................................."
    !write(*,*) "connected_to_start =", this%connected_to_start

    tcfc = COUNT(connected,1) ! count the total number of TRUE (whole domain)
    allocate(this%unique_connected_cell_ids_glb(tcfc))
    ! how do I check if a glb id is owned by my process and what the corresponding
    ! local id is?
    j = 0
    do k = 1,res_grid%nlmax
      glb_id = res_grid%nG2A(res_grid%nL2G(k))
      ! count how many of the local ids I own have a gbl id that is connected
      do i=1,size(this%allfrac_unique_cell_ids_glb)
        if (this%allfrac_unique_cell_ids_glb(i) == glb_id) exit
      enddo
      if (connected(i)) j = j+1
    enddo
    allocate(this%unique_connected_cell_ids(j)) ! ok if j=0
    if (size(this%unique_connected_cell_ids) > 0) then
      j = 0
      do k = 1,res_grid%nlmax
        glb_id = res_grid%nG2A(res_grid%nL2G(k))
        ! load the local ids I own which have a gbl id that is connected
        do i=1,size(this%allfrac_unique_cell_ids_glb)
          if (this%allfrac_unique_cell_ids_glb(i) == glb_id) exit
        enddo
        if (connected(i)) then
          j = j+1
          this%unique_connected_cell_ids(j) = k ! k is the local id
        endif
      enddo
    endif

    ! if in parallel, the list in this%unique_connected_cell_ids are local to
    ! the processor, so we need to make a global id list from all processors
    allocate(uniq_conn_num_cells(option%comm%size))
    j = size(this%unique_connected_cell_ids)
    call MPI_Allgather(j,ONE_INTEGER,MPI_INTEGER,uniq_conn_num_cells, &
                       ONE_INTEGER,MPI_INTEGER,option%mycomm,ierr);CHKERRQ(ierr)
    call MPI_Barrier(option%mycomm,ierr);CHKERRQ(ierr)
    allocate(temp_gbl_cell_ids(j)) ! ok if j=0
    allocate(this%unique_connected_cell_ids_glb(SUM(uniq_conn_num_cells)))
    do i=1,j ! never enters if j=0
      temp_gbl_cell_ids(i) = &
                res_grid%nG2A(res_grid%nL2G(this%unique_connected_cell_ids(i)))
    enddo
    do i = 1,option%comm%size
      if (i == 1) then
        cell_positions(i) = 0
      else
        cell_positions(i) = cell_positions(i-1) + uniq_conn_num_cells(i-1)
      endif
    enddo
    if (SUM(uniq_conn_num_cells) > 0) then
      call MPI_Gatherv(temp_gbl_cell_ids,j,MPI_INTEGER, &
                       this%unique_connected_cell_ids_glb,uniq_conn_num_cells, &
                       cell_positions,MPI_INTEGER,ZERO_INTEGER, &
                       option%mycomm,ierr);CHKERRQ(ierr)
      ! only rank zero knows this%unique_connected_cell_ids_glb, so share it
      call MPI_Barrier(option%mycomm,ierr);CHKERRQ(ierr)
      call MPI_Bcast(this%unique_connected_cell_ids_glb, &
                     SUM(uniq_conn_num_cells),MPI_INTEGER,ZERO_INTEGER, &
                     option%mycomm,ierr);CHKERRQ(ierr)
    endif
    deallocate(uniq_conn_num_cells)
    deallocate(temp_gbl_cell_ids)

    !write(*,*) "............................................................."
    !write(*,*) "connected cell ids =", this%unique_connected_cell_ids_glb

    if (tcfc < size(this%allfrac_unique_cell_ids_glb)) then
      option%io_buffer = 'GEOTHERMAL_FRACTURE_MODEL: There are disconnected &
        &fractures in the domain. Disconnected fractures have been clipped.'
      call PrintMsg(option)
    endif

    ! After clipping, the number of fracture intersections may be reduced.
    ! Need to keep only the cells inside this%frac_intersection_ids that also
    ! still exist in this%unique_connected_cell_ids_glb. Reuse temp_gbl_cell_ids.
    allocate(temp_gbl_cell_ids(size(this%frac_intersection_ids)))
    temp_gbl_cell_ids(:) = UNINITIALIZED_INTEGER
    k = 0
    do i = 1, size(this%frac_intersection_ids)
      if (any(this%unique_connected_cell_ids_glb == &
              this%frac_intersection_ids(i))) then
        k = k+1
        temp_gbl_cell_ids(k) = this%frac_intersection_ids(i)
      endif
    enddo
    deallocate(this%frac_intersection_ids)
    allocate(this%frac_intersection_ids(k)) ! resize (shorten)
    this%frac_intersection_ids(:) = UNINITIALIZED_INTEGER
    if (k > 0) this%frac_intersection_ids = temp_gbl_cell_ids(1:k)
    deallocate(temp_gbl_cell_ids)
  endif ! if (this%clip_unconnected_fracs) then

  deallocate(cell_positions)
  ! filename needs to be consistant w/ filename located at end of PMFracPrintStats()
  filename = trim(this%option%global_prefix)//'.info' ! ejb: can call our other file .fracinfo
  fid = 19181723    ! ejb: can randomly choose an fid
110 format(100es16.8) ! ejb: 110 <-- shortcut, wherever you see "110", insert the formatting code
  open(unit=fid,file=filename,action="write",status="replace") ! ejb: open and write are fortran inherent fxns, best to all-caps them if possible
  write(fid,'(a)',advance="yes") "FRACTURE NETWORK METRICS FOR PFLOTRAN &
                        &INPUT DECK: "//trim(this%option%global_prefix)//".in"
  call GETCWD(path,status) ! ejb: another fortran-inherent
  if (status==0) write(fid,'(a)',advance="yes") "INPUT DECK PATH: " &
                                                //trim(path)
  if (status/=0) write(fid,'(a)',advance="yes") "INPUT DECK PATH: unknown"
  write(fid,'(a)',advance="yes") "TIMESTAMP OF FILE CREATION: "//FDATE()
  write(fid,'(a)',advance="yes") " "
  write(fid,'(a)',advance="no") "NUMBER OF GRID CELLS IN DOMAIN: "
  write(fid,'(I10)',advance="yes") res_grid%nmax
  
  if (this%clip_unconnected_fracs) then
    write(fid,'(a)',advance="no") "NUMBER OF FRACTURED GRID CELLS: "
    write(fid,'(I10)',advance="yes") size(this%unique_connected_cell_ids_glb)
    write(fid,'(a)',advance="no") "PERCENT [%] OF DOMAIN GRID CELLS THAT &
                                   &CONTAIN FRACTURES: "
    write(fid,110,advance="yes") &
                         (REAL(size(this%unique_connected_cell_ids_glb))/ &
                          REAL(res_grid%nmax))*100.d0
    L22 = 0.d0
    do k = 1,size(this%unique_connected_cell_ids)
      material_auxvar => this%realization%patch%aux%material% &
                    auxvars(res_grid%nL2G(this%unique_connected_cell_ids(k)))
      L22 = L22 + (2.d0*((material_auxvar%volume)**(1.d0/3.d0))**(2)) ! [m^2]
    enddo
  else
    write(fid,'(a)',advance="no") "NUMBER OF FRACTURED GRID CELLS: "
    write(fid,'(I10)',advance="yes") size(this%allfrac_unique_cell_ids_glb)
    write(fid,'(a)',advance="no") "PERCENT [%] OF DOMAIN GRID CELLS THAT &
                                   &CONTAIN FRACTURES: "
    write(fid,110,advance="yes") &
      (REAL(size(this%allfrac_unique_cell_ids_glb))/REAL(res_grid%nmax))*100.d0
    L22 = 0.d0
    do k = 1,size(this%allfrac_unique_cell_ids)
      material_auxvar => this%realization%patch%aux%material% &
                  auxvars(res_grid%nL2G(this%allfrac_unique_cell_ids(k)))
      L22 = L22 + (2.d0*((material_auxvar%volume)**(1.d0/3.d0))**(2)) ! [m^2]
    enddo
  endif

  ! need to get all processes value for L22 and sum them
  call MPI_Allreduce(MPI_IN_PLACE,L22,option%comm%size,MPI_DOUBLE_PRECISION, &
                     MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
  write(fid,'(a)',advance="no") "SURFACE AREA [m^2] OF FRACTURES: "
  write(fid,110,advance="yes") L22
  write(fid,'(a)',advance="no") "NUMBER OF FRAC INTERSECTIONS IN DOMAIN: "
  write(fid,'(I10)',advance="yes") size(this%frac_intersection_ids)

  write(fid,'(a)',advance="yes") ""
  ! header for time varying fracture network characteristics
  write(fid,'(a)',advance="no") "TIME ["//trim(this%output_option%tunit)//"],"
  write(fid,'(a)',advance="no") "VOL [m^3] OF FRACTURES,"
  write(fid,'(a)',advance="no") "MIN FRAC APERTURE [m],"
  write(fid,'(a)',advance="no") "MAX FRAC APERTURE [m],"
  write(fid,'(a)',advance="no") "RANGE FRAC APERTURE [m],"
  write(fid,'(a)',advance="no") "MIDRANGE FRAC APERTURE [m],"
  write(fid,'(a)',advance="no") "MEAN FRAC APERTURE [m],"
  write(fid,'(a)',advance="no") "VARIANCE FRAC APERTURE [m],"
  write(fid,'(a)',advance="no") "POP SKEWNESS FRAC APERTURE [m],"
  write(fid,'(a)',advance="no") "SAMPLE SKEWNESS FRAC APERTURE [m],"
  write(fid,'(a)',advance="no") "MIN FRAC PERM-X [m],"
  write(fid,'(a)',advance="no") "MAX FRAC PERM-X [m],"
  write(fid,'(a)',advance="no") "RANGE FRAC PERM-X [m],"
  write(fid,'(a)',advance="no") "MIDRANGE FRAC PERM-X [m],"
  write(fid,'(a)',advance="no") "MEAN FRAC PERM-X [m],"
  write(fid,'(a)',advance="no") "VARIANCE FRAC PERM-X [m],"
  write(fid,'(a)',advance="no") "POP SKEWNESS FRAC PERM-X [m],"
  write(fid,'(a)',advance="no") "SAMPLE SKEWNESS FRAC PERM-X [m],"
  write(fid,'(a)',advance="no") "MIN FRAC PERM-Y [m],"
  write(fid,'(a)',advance="no") "MAX FRAC PERM-Y [m],"
  write(fid,'(a)',advance="no") "RANGE FRAC PERM-Y [m],"
  write(fid,'(a)',advance="no") "MIDRANGE FRAC PERM-Y [m],"
  write(fid,'(a)',advance="no") "MEAN FRAC PERM-Y [m],"
  write(fid,'(a)',advance="no") "VARIANCE FRAC PERM-Y [m],"
  write(fid,'(a)',advance="no") "POP SKEWNESS FRAC PERM-Y [m],"
  write(fid,'(a)',advance="no") "SAMPLE SKEWNESS FRAC PERM-Y [m],"
  write(fid,'(a)',advance="no") "MIN FRAC PERM-Z [m],"
  write(fid,'(a)',advance="no") "MAX FRAC PERM-Z [m],"
  write(fid,'(a)',advance="no") "RANGE FRAC PERM-Z [m],"
  write(fid,'(a)',advance="no") "MIDRANGE FRAC PERM-Z [m],"
  write(fid,'(a)',advance="no") "MEAN FRAC PERM-Z [m],"
  write(fid,'(a)',advance="no") "VARIANCE FRAC PERM-Z [m],"
  write(fid,'(a)',advance="no") "POP SKEWNESS FRAC PERM-Z [m],"
  write(fid,'(a)',advance="no") "SAMPLE SKEWNESS FRAC PERM-Z [m],"
  close(fid)

  ! filename_frac_info needs to be consistent w/ filename_frac_info located at end of PMFracPrintStats()
  ! EJB: filename for printing out fracture geometry information: 
  filename_frac_info = trim(this%option%global_prefix)//'.frac_geom_info'
  fid_frac_info = 20846713
  OPEN(unit=fid_frac_info,file=filename_frac_info,action="write",status="replace") 
  WRITE(fid_frac_info,'(a)',advance="yes") "FRACTURE NETWORK METRICS FOR PFLOTRAN &
                        &INPUT DECK: "//trim(this%option%global_prefix)//".in"
  call GETCWD(path,status) ! ejb: another fortran-inherent
  if (status==0) WRITE(fid_frac_info,'(a)',advance="yes") "INPUT DECK PATH: " &
                                                //trim(path)
  if (status/=0) WRITE(fid_frac_info,'(a)',advance="yes") "INPUT DECK PATH: unknown"
  WRITE(fid_frac_info,'(a)',advance="yes") "TIMESTAMP OF FILE CREATION: "//FDATE()
  WRITE(fid_frac_info,'(a)',advance="no") "fracture_id,"
  WRITE(fid_frac_info,'(a)',advance="no") "center_x,center_y,center_z,"
  WRITE(fid_frac_info,'(a)',advance="no") "normal_x,normal_y,normal_z,"
  WRITE(fid_frac_info,'(a)',advance="no") "radius_x,radius_y,radius_z,"
  WRITE(fid_frac_info,'(a)',advance="no") "aperture_initial,"
  WRITE(fid_frac_info,'(a)',advance="no") "strike_x,strike_y,strike_z,"
  WRITE(fid_frac_info,'(a)',advance="yes") "dip"
  ! WRITE(fid_frac_info,'(a)',advance="yes") ""
  CLOSE(fid_frac_info)

  ! now, call our subroutine to write frac geom stats: 
  call PMFracGeomPrintStats(this)

  if (.not.associated(this%realization%reaction)) then
    this%update_material_auxvar_perm = PETSC_TRUE
  else
    if (.not.this%realization%reaction%update_permeability) then
      this%update_material_auxvar_perm = PETSC_TRUE
    endif
  endif

  allocate(this%max_frac_kx(res_grid%nmax))
  allocate(this%max_frac_ky(res_grid%nmax))
  allocate(this%max_frac_kz(res_grid%nmax))
  this%max_frac_kx = 0.d0; this%max_frac_ky = 0.d0; this%max_frac_kz = 0.d0
  allocate(this%max_frac_hap(res_grid%nmax))
  this%max_frac_hap = 0.d0

  allocate(this%sum_frac_kx(res_grid%nmax))
  allocate(this%sum_frac_ky(res_grid%nmax))
  allocate(this%sum_frac_kz(res_grid%nmax))
  this%sum_frac_kx = 0.d0; this%sum_frac_ky = 0.d0; this%sum_frac_kz = 0.d0
  allocate(this%sum_frac_hap(res_grid%nmax))
  this%sum_frac_hap = 0.d0

  deallocate(temp_cell_ids)
  deallocate(temp_allfrac_cell_ids)
  deallocate(d_vertex)

end subroutine PMFracSetup

! ************************************************************************** !

subroutine PMFracClipCells(this,connected_to_region,total_sum_glb, &
                           cell_positions,region_name)
  !
  ! Figures out what cells are in the clip region and connected to the region input
  !
  use Grid_module
  use Region_module

  implicit none

  class(pm_fracture_type) :: this
  PetscBool, pointer :: connected_to_region(:)
  PetscInt :: total_sum_glb
  PetscInt, pointer :: cell_positions(:)
  character(len=MAXWORDLENGTH) :: region_name

  type(option_type), pointer :: option
  type(region_type), pointer :: region
  type(grid_type), pointer :: res_grid
  PetscInt :: k,i
  PetscInt :: myregion_num_cells,region_num_cells
  PetscInt, pointer :: allregion_num_cells(:)
  PetscInt, pointer :: myregion_cell_ids(:)
  PetscBool, pointer :: connected(:)
  PetscErrorCode :: ierr
  PetscInt, pointer :: visited_cells(:)

  option => this%option
  res_grid => this%realization%patch%grid

  region => RegionGetPtrFromList(region_name, &
                                 this%realization%patch%region_list)

  if (.not.associated(region)) then
    option%io_buffer = 'REGION ' // trim(region_name) // ' in &
                        &CLIP_UNCONNECTED_FRACTURES,REGION not found among the &
                        &defined REGIONs in the input deck.'
    call PrintErrMsg(option)
  endif

  allocate(allregion_num_cells(option%comm%size))
  myregion_num_cells = region%num_cells
  if (myregion_num_cells > 0) then
    allocate(myregion_cell_ids(myregion_num_cells))
    do i = 1,myregion_num_cells
      ! get global id for the region cell ids
      myregion_cell_ids(i) = &
        res_grid%nG2A(res_grid%nL2G(region%cell_ids(i)))
    enddo
  else
    allocate(myregion_cell_ids(1))
    myregion_cell_ids = UNINITIALIZED_INTEGER
  endif

  call MPI_Barrier(option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Allgather(myregion_num_cells,ONE_INTEGER,MPI_INTEGER, &
                     allregion_num_cells,ONE_INTEGER,MPI_INTEGER, &
                     option%mycomm,ierr);CHKERRQ(ierr)

  region_num_cells = SUM(allregion_num_cells)
  do i = 1,option%comm%size
    if (i == 1) then
      cell_positions(i) = 0
    else
      cell_positions(i) = cell_positions(i-1) + allregion_num_cells(i-1)
    endif
  enddo

  allocate(this%clip_region_cell_ids(region_num_cells))
  call MPI_Gatherv(myregion_cell_ids,myregion_num_cells,MPI_INTEGER, &
                   this%clip_region_cell_ids,allregion_num_cells, &
                   cell_positions,MPI_INTEGER,ZERO_INTEGER, &
                   option%mycomm,ierr);CHKERRQ(ierr)
  ! only rank zero knows this%clip_region_cell_ids, so we need to share it
  call MPI_Barrier(option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Bcast(this%clip_region_cell_ids,region_num_cells,MPI_INTEGER, &
                 ZERO_INTEGER,option%mycomm,ierr);CHKERRQ(ierr)
  deallocate(myregion_cell_ids)
  deallocate(allregion_num_cells)
  nullify(connected)

  allocate(connected_to_region(total_sum_glb))
  connected_to_region(:) = PETSC_FALSE
  allocate(visited_cells(total_sum_glb))

  do k = 1,size(this%clip_region_cell_ids)
    ! this is a recursive subroutine which will call itself
    visited_cells = 0
    call PMFracConnectedness(this,this%clip_region_cell_ids(k),visited_cells, &
                             connected_to_region)
  enddo

end subroutine PMFracClipCells

! ************************************************************************** !

recursive subroutine PMFracConnectedness(this,natural_id,visited_cells, &
                                         connected)
  !
  ! Figures out the fracture cells which are connected between a start region
  ! and an end region.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 04/23/2025
  !
  use Grid_module
  use Connection_module

  implicit none

  class(pm_fracture_type) :: this
  PetscInt :: natural_id
  PetscInt, pointer :: visited_cells(:)

  type(grid_type), pointer :: grid
  type(connection_set_type), pointer :: cur_connection_set
  type(connection_set_list_type), pointer :: connection_set_list
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: natural_id_up, natural_id_dn
  PetscInt :: iconn, i
  PetscBool, pointer :: connected(:)

  grid => this%realization%patch%grid
  connection_set_list => grid%internal_connection_set_list

  if (.not.this%clip_unconnected_fracs) return








  ! something with the connection set isn't working. connection sets are local
  ! so how do i make the routine work to check for connections across processors?









  if (any(this%allfrac_unique_cell_ids_glb == natural_id)) then
    ! note: FINDLOC only works with gfortran > 9.x
    !i = FINDLOC(this%allfrac_unique_cell_ids_glb,natural_id,1)
    do i=1,size(this%allfrac_unique_cell_ids_glb)
      if (this%allfrac_unique_cell_ids_glb(i) == natural_id) exit
   enddo
   if (visited_cells(i) == 1) then
      return
   else
      visited_cells(i) = 1
   endif
   connected(i) = PETSC_TRUE
  else
    return
  endif

  cur_connection_set => connection_set_list%first
  do
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn)
      natural_id_up = grid%nG2A(ghosted_id_up)
      natural_id_dn = grid%nG2A(ghosted_id_dn)

      if (natural_id_up == natural_id) then
        !write(*,*) " "
        !write(*,*) "Cell",natural_id,"is the natural_id_up of connection", iconn
        !write(*,*) "Connection", iconn, " connects cell",natural_id,"with cell", &
        !             natural_id_dn, "."
        if (any(this%allfrac_unique_cell_ids_glb == natural_id_dn)) then
          !write(*,*) "!! Cell", natural_id_dn, "is a fracture cell!!"
          ! note: FINDLOC only works with gfortran > 9.x
          !i = FINDLOC(this%allfrac_unique_cell_ids_glb,natural_id_dn,1)
          do i=1,size(this%allfrac_unique_cell_ids_glb)
            if (this%allfrac_unique_cell_ids_glb(i) == natural_id) exit
          enddo
          if (connected(i)) then
            call PMFracConnectedness(this,natural_id_dn,visited_cells,connected)
          endif
        endif
      endif

      if (natural_id_dn == natural_id) then
        !write(*,*) " "
        !write(*,*) "Cell",natural_id,"is the natural_id_dn of connection", iconn
        !write(*,*) "Connection", iconn, " connects cell",natural_id,"with cell", &
        !             natural_id_up, "."
        if (any(this%allfrac_unique_cell_ids_glb == natural_id_up)) then
          !write(*,*) "!! Cell", natural_id_up, "is a fracture cell!!"
          ! note: FINDLOC only works with gfortran > 9.x
          !i = FINDLOC(this%allfrac_unique_cell_ids_glb,natural_id_up,1)
          do i=1,size(this%allfrac_unique_cell_ids_glb)
            if (this%allfrac_unique_cell_ids_glb(i) == natural_id) exit
          enddo

          if (connected(i)) then
            call PMFracConnectedness(this,natural_id_up,visited_cells,connected)
          endif
        endif
      endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo

end subroutine PMFracConnectedness

! ************************************************************************** !

subroutine PMFracSetRealization(this)
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !
  use Realization_Subsurface_class

  implicit none

  class(pm_fracture_type) :: this

  this%realization => RealizationCast(this%realization_base)

end subroutine PMFracSetRealization

! ************************************************************************** !

subroutine PMFracInitializeRun(this)
  !
  ! Initializes the run associated with the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !
  use Material_Aux_module
  use TH_Aux_module
  use General_Aux_module
  use Output_Aux_module
  use Field_module
  use Grid_module

  implicit none

  class(pm_fracture_type) :: this

  type(material_auxvar_type), pointer :: material_auxvar
  type(th_auxvar_type), pointer :: th_auxvar
  type(general_auxvar_type), pointer :: general_auxvar
  type(output_variable_list_type), pointer :: output_var_list
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(fracture_type), pointer :: cur_fracture
  PetscReal, pointer :: perm0_xx_p(:), perm0_yy_p(:), perm0_zz_p(:)
  PetscInt :: icell,k
  PetscReal :: kx,ky,kz,L
  PetscBool :: sit1,sit2
  PetscInt, pointer :: cellid_list(:)
  PetscErrorCode :: ierr

  option => this%option
  field => this%realization%field
  grid => this%realization%patch%grid

  call VecGetArrayRead(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)
  this%max_frac_kx = 0.d0; this%max_frac_ky = 0.d0; this%max_frac_kz = 0.d0
  this%sum_frac_kx = 0.d0; this%sum_frac_ky = 0.d0; this%sum_frac_kz = 0.d0
  this%max_frac_hap = 0.d0; this%sum_frac_hap = 0.d0

  sit1 = .not.this%clip_unconnected_fracs
  sit2 = PETSC_FALSE

  cur_fracture => this%fracture_list
  do
    if (.not.associated(cur_fracture)) exit

    do k = 1,cur_fracture%ncells
      icell = cur_fracture%cell_ids(k) ! local cell ids
      material_auxvar => &
        this%realization%patch%aux%material%auxvars(grid%nL2G(icell))
      ! if cells aren't too pancaked, L should be a good estimate of length
      L = (material_auxvar%volume)**(1.d0/3.d0)

      cur_fracture%kx0(k) = perm0_xx_p(icell)
      cur_fracture%ky0(k) = perm0_yy_p(icell)
      cur_fracture%kz0(k) = perm0_zz_p(icell)

      call PMFracCalcK(L,cur_fracture%hap0,cur_fracture%RM,kx,ky,kz)

      ! record the maximum permeability encountered so far
      this%max_frac_kx(icell) = max(kx,this%max_frac_kx(icell))
      this%max_frac_ky(icell) = max(ky,this%max_frac_ky(icell))
      this%max_frac_kz(icell) = max(kz,this%max_frac_kz(icell))
      this%max_frac_hap(icell) = &
        max(cur_fracture%hap0,this%max_frac_hap(icell))

      ! record the sum of permeability encountered so far
      this%sum_frac_kx(icell) = this%sum_frac_kx(icell) + kx
      this%sum_frac_ky(icell) = this%sum_frac_ky(icell) + ky
      this%sum_frac_kz(icell) = this%sum_frac_kz(icell) + kz
      this%sum_frac_hap(icell) = this%sum_frac_hap(icell) + cur_fracture%hap0

      if (this%clip_unconnected_fracs) then
        sit2 = any(icell == this%unique_connected_cell_ids)
      endif

      if (sit1 .or. sit2) then
        material_auxvar%porosity_0 = material_auxvar%porosity_0 + &
                                     (cur_fracture%hap0/L)
        material_auxvar%porosity_base = material_auxvar%porosity_base + &
                                        (cur_fracture%hap0/L)
        material_auxvar%porosity = material_auxvar%porosity + &
                                   (cur_fracture%hap0/L)
      endif

      if (option%iflowmode == G_MODE) then
        general_auxvar => &
          this%realization%patch%aux%General%auxvars(1,grid%nL2G(icell))
        cur_fracture%prev_temperature(k) = general_auxvar%temp
      elseif (option%iflowmode == TH_MODE) then
        th_auxvar => &
          this%realization%patch%aux%TH%auxvars(grid%nL2G(icell))
        cur_fracture%prev_temperature(k) = th_auxvar%temp
      endif
    enddo

    cur_fracture => cur_fracture%next
  enddo

  if (this%clip_unconnected_fracs) then
    cellid_list => this%unique_connected_cell_ids ! local cell ids
  else
    cellid_list => this%allfrac_unique_cell_ids ! local cell ids
  endif

  ! update domain material permeability of cells with fractures
  if (associated(cellid_list)) then
    do k = 1,size(cellid_list)
      icell = cellid_list(k)
      material_auxvar => &
        this%realization%patch%aux%material%auxvars(grid%nL2G(icell))
      if (this%frac_intersection_type == 2) then ! maximum
        perm0_xx_p(icell) = perm0_xx_p(icell) + this%max_frac_kx(icell)
        perm0_yy_p(icell) = perm0_yy_p(icell) + this%max_frac_ky(icell)
        perm0_zz_p(icell) = perm0_zz_p(icell) + this%max_frac_kz(icell)
        material_auxvar%fracture_aperture = this%max_frac_hap(icell)
      endif
      if (this%frac_intersection_type == 1) then ! summation
        perm0_xx_p(icell) = perm0_xx_p(icell) + this%sum_frac_kx(icell)
        perm0_yy_p(icell) = perm0_yy_p(icell) + this%sum_frac_ky(icell)
        perm0_zz_p(icell) = perm0_zz_p(icell) + this%sum_frac_kz(icell)
        material_auxvar%fracture_aperture = this%sum_frac_hap(icell)
      endif
    enddo
  endif

  if (this%update_material_auxvar_perm) then
    if (associated(cellid_list)) then
      do k = 1,size(cellid_list)
        icell = cellid_list(k)
        material_auxvar => &
          this%realization%patch%aux%material%auxvars(grid%nL2G(icell))
        material_auxvar%permeability(1) = perm0_xx_p(icell)
        material_auxvar%permeability(2) = perm0_yy_p(icell)
        material_auxvar%permeability(3) = perm0_zz_p(icell)
      enddo
    endif
  endif

  call VecRestoreArrayRead(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)

  ! Setup the output variables for snapshot and observation files
  output_var_list => this%realization%output_option%output_snap_variable_list
  call PMFracSetPlotVariables(output_var_list,this)
  if (.not.associated( &
                this%realization%output_option%output_snap_variable_list, &
                this%realization%output_option%output_obs_variable_list)) then
    output_var_list => this%realization%output_option%output_obs_variable_list
    call PMFracSetPlotVariables(output_var_list,this)
  endif

end subroutine PMFracInitializeRun

! ************************************************************************** !

subroutine PMFracFinalizeRun(this)
  !
  ! Finalizes the run associated with the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !

  implicit none

  class(pm_fracture_type) :: this

end subroutine PMFracFinalizeRun

! ************************************************************************** !

subroutine PMFracInitializeTimestep(this)
  !
  ! Initializes the timestep associated with the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !
  use Material_Aux_module
  use Field_module
  use Grid_module

  implicit none

  class(pm_fracture_type) :: this

  type(material_auxvar_type), pointer :: material_auxvar
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(fracture_type), pointer :: cur_fracture
  PetscReal, pointer :: perm0_xx_p(:), perm0_yy_p(:), perm0_zz_p(:)
  PetscInt :: icell,k
  PetscErrorCode :: ierr

  field => this%realization%field
  grid => this%realization%patch%grid

  call VecGetArrayRead(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)

  cur_fracture => this%fracture_list
  do
    if (.not.associated(cur_fracture)) exit

    do k = 1,cur_fracture%ncells
      icell = cur_fracture%cell_ids(k)
      material_auxvar => &
        this%realization%patch%aux%material%auxvars(grid%nL2G(icell))
      ! reset permeability back to initial permeability
      perm0_xx_p(icell) = cur_fracture%kx0(k)
      perm0_yy_p(icell) = cur_fracture%ky0(k)
      perm0_zz_p(icell) = cur_fracture%kz0(k)
    enddo

    cur_fracture => cur_fracture%next
  enddo

  call VecRestoreArrayRead(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)

end subroutine PMFracInitializeTimestep

! ************************************************************************** !

subroutine PMFracFinalizeTimestep(this)
  !
  ! Finalizes the timestep associated with the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !

  implicit none

  class(pm_fracture_type) :: this

end subroutine PMFracFinalizeTimestep

! ************************************************************************** !

subroutine PMFracReadPMBlock(this,input)
  !
  ! Reads input file parameters associated with the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023

  implicit none

  class(pm_fracture_type) :: this
  type(input_type), pointer :: input

  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found,fracfam_given,frac_given

  option => this%option
  fracfam_given = PETSC_FALSE
  frac_given = PETSC_FALSE
  input%ierr = INPUT_ERROR_NONE
  error_string = 'GEOTHERMAL_FRACTURE_MODEL'

  option%io_buffer = 'pflotran card:: GEOTHERMAL_FRACTURE_MODEL'
  call PrintMsg(option)

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    found = PETSC_FALSE

    ! Read keywords within GEOTHERMAL_FRACTURE_MODEL block:
    select case(trim(word))
    !-------------------------------------
      case('THERMAL_EXPANSION_COEFFICIENT')
        call InputReadDouble(input,option,this%t_coeff)
        call InputErrorMsg(input,option,'THERMAL_EXPANSION_COEFFICIENT', &
                         error_string)
        cycle
    !-------------------------------------
      case('MAXIMUM_NUMBER_OF_FRACTURES')
        call InputReadInt(input,option,this%max_frac)
        call InputErrorMsg(input,option,'MAXIMUM_NUMBER_OF_FRACTURES', &
                         error_string)
        cycle
    !-------------------------------------
      case('FRACTURE_INTERSECTION_PERM')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,word,error_string)
        call StringToUpper(word)
        select case (trim(word))
          case('MAX','MAXIMUM')
            this%frac_intersection_type = 2
          case('SUM','SUMMATION')
            this%frac_intersection_type = 1
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        end select
        cycle
    !-------------------------------------

    end select

    ! Read sub-blocks within GEOTHERMAL_FRACTURE_MODEL block:
    error_string = 'GEOTHERMAL_FRACTURE_MODEL'
    call PMFracReadFracture(this,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'GEOTHERMAL_FRACTURE_MODEL'
    call PMFracReadFracFam(this,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'GEOTHERMAL_FRACTURE_MODEL'
    call PMFracReadClipRegion(this,input,option,word,error_string,found)
    if (found) cycle

    if (.not. found) then
      option%io_buffer = 'Keyword "' // trim(word) // &
                         '" does not exist for GEOTHERMAL_FRACTURE_MODEL.'
      call PrintErrMsg(option)
    endif

  enddo
  call InputPopBlock(input,option)

if (Uninitialized(this%t_coeff)) then
  option%io_buffer = 'THERMAL_EXPANSION_COEFFICIENT not provided in &
    &the GEOTHERMAL_FRACTURE_MODEL block.'
  call PrintErrMsg(option)
endif

if (Uninitialized(this%max_frac)) then
  this%max_frac = 100000000  ! set to a very high number so the max
                             ! is never reached
endif

if (associated(this%fracfam_list)) fracfam_given = PETSC_TRUE
if (associated(this%fracture_list)) frac_given = PETSC_TRUE

if (.not.fracfam_given .and. .not.frac_given) then
  option%io_buffer = 'At least one FRACTURE or FRACTURE_FAMILY block must &
    &be given in the GEOTHERMAL_FRACTURE_MODEL block.'
  call PrintErrMsg(option)
endif

end subroutine PMFracReadPMBlock

! ************************************************************************** !

subroutine PMFracReadFracture(pm_fracture,input,option,keyword,error_string, &
                              found)
  !
  ! Reads input file parameters associated with the fracture model grid.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !

  implicit none

  class(pm_fracture_type) :: pm_fracture
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  character(len=MAXWORDLENGTH) :: word
  type(fracture_type), pointer :: new_fracture,cur_fracture
  PetscReal :: vmag
  PetscInt :: num_errors
  PetscBool :: added

  error_string = trim(error_string) // ',FRACTURE'
  found = PETSC_TRUE
  added = PETSC_FALSE
  num_errors = 0

  select case(trim(keyword))
  !-------------------------------------
    case('FRACTURE')
      allocate(new_fracture)
      new_fracture => FractureCreate()
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case(trim(word))
        !-----------------------------
          case('ID')
            call InputReadInt(input,option,new_fracture%id)
            call InputErrorMsg(input,option,'ID',error_string)
        !-----------------------------
          case('HYDRAULIC_APERTURE')
            call InputReadDouble(input,option,new_fracture%hap0)
            call InputErrorMsg(input,option,'HYDRAULIC_APERTURE',error_string)
        !-----------------------------
          case('RADIUS_X')
            call InputReadDouble(input,option,new_fracture%radx)
            call InputErrorMsg(input,option,'RADIUS_X',error_string)
        !-----------------------------
          case('RADIUS_Y')
            call InputReadDouble(input,option,new_fracture%rady)
            call InputErrorMsg(input,option,'RADIUS_Y',error_string)
        !-----------------------------
          case('RADIUS_Z')
            call InputReadDouble(input,option,new_fracture%radz)
            call InputErrorMsg(input,option,'RADIUS_Z',error_string)
        !-----------------------------
          case('CENTER')
            call InputReadDouble(input,option,new_fracture%center%x)
            call InputErrorMsg(input,option,'CENTER, X COORDINATE', &
                             error_string)
            call InputReadDouble(input,option,new_fracture%center%y)
            call InputErrorMsg(input,option,'CENTER, Y COORDINATE', &
                             error_string)
            call InputReadDouble(input,option,new_fracture%center%z)
            call InputErrorMsg(input,option,'CENTER, Z COORDINATE', &
                             error_string)
        !-----------------------------
          case('NORMAL_VECTOR')
            call InputReadDouble(input,option,new_fracture%normal%x)
            call InputErrorMsg(input,option,'NORMAL_VECTOR, X COORDINATE', &
                             error_string)
            call InputReadDouble(input,option,new_fracture%normal%y)
            call InputErrorMsg(input,option,'NORMAL_VECTOR, Y COORDINATE', &
                             error_string)
            call InputReadDouble(input,option,new_fracture%normal%z)
            call InputErrorMsg(input,option,'NORMAL_VECTOR, Z COORDINATE', &
                             error_string)
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)
      !------ error messaging ----------------------------------------
      if (Uninitialized(new_fracture%id)) then
        option%io_buffer = 'ERROR: ID must be specified in ' // &
                           trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracture%hap0)) then
        option%io_buffer = 'ERROR: HYDRAULIC_APERTURE must be specified in ' // &
                           trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracture%radx)) then
        option%io_buffer = 'ERROR: RADIUS_X must be specified in ' // &
                           trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracture%rady)) then
        option%io_buffer = 'ERROR: RADIUS_Y must be specified in ' // &
                           trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracture%radz)) then
        option%io_buffer = 'ERROR: RADIUS_Z must be specified in ' // &
                           trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracture%center%x) .or. &
          Uninitialized(new_fracture%center%y) .or. &
          Uninitialized(new_fracture%center%z)) then
        option%io_buffer = 'ERROR: CENTER coordinate must be specified in ' // &
                           trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracture%normal%x) .or. &
          Uninitialized(new_fracture%normal%y) .or. &
          Uninitialized(new_fracture%normal%z)) then
        option%io_buffer = 'ERROR: NORMAL_VECTOR coordinate must be specified &
                            &in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      else
        vmag = sqrt((new_fracture%normal%x)**2.d0 + &
                  (new_fracture%normal%y)**2.d0 + &
                  (new_fracture%normal%z)**2.d0)
        new_fracture%normal%x = new_fracture%normal%x/vmag
        new_fracture%normal%y = new_fracture%normal%y/vmag
        new_fracture%normal%z = new_fracture%normal%z/vmag
      endif
      if (num_errors > 0) then
        write(option%io_buffer,*) num_errors
        option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                           &the FRACTURE block. See above to fix.'
        call PrintErrMsg(option)
      endif
      !------ add fracture to list -----------------------------------
      if (.not.associated(pm_fracture%fracture_list)) then
        pm_fracture%fracture_list => new_fracture
      else
        cur_fracture => pm_fracture%fracture_list
        do
          if (.not.associated(cur_fracture)) exit
          if (.not.associated(cur_fracture%next)) then
            cur_fracture%next => new_fracture
            added = PETSC_TRUE
          endif
          if (added) exit
          cur_fracture => cur_fracture%next
        enddo
      endif
      nullify(new_fracture)
  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

end subroutine PMFracReadFracture

! ************************************************************************** !

subroutine PMFracReadFracFam(pm_fracture,input,option,keyword,error_string, &
                            found)
  !
  ! Reads input file parameters associated with the fracture model grid.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/19/2024
  !

  implicit none

  class(pm_fracture_type) :: pm_fracture
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  type(fracfam_type), pointer :: new_fracfam, cur_fracfam
  character(len=MAXWORDLENGTH) :: word
  PetscBool :: added
  PetscInt :: num_errors

  error_string = trim(error_string) // ',FRACTURE_FAMILY'
  found = PETSC_TRUE
  added = PETSC_FALSE
  num_errors = 0

  select case(trim(keyword))
  !-------------------------------------
    case('FRACTURE_FAMILY')
      allocate(new_fracfam)
      new_fracfam => FractureFamCreate()
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case(trim(word))
        !-----------------------------
          case('ID')
            call InputReadInt(input,option,new_fracfam%id)
            call InputErrorMsg(input,option,'ID',error_string)
        !-----------------------------
          case('NUMBER_OF_FRACTURES')
            call InputReadInt(input,option,new_fracfam%intensity_fracfam)
            call InputErrorMsg(input,option,'NUMBER_OF_FRACTURES',error_string)
        !-----------------------------
          case('HYDRAULIC_APERTURE')
            call InputPushBlock(input,option)
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                call InputReadCard(input,option,word)
                call InputErrorMsg(input,option,'keyword',error_string)
                call StringToUpper(word)
                select case(trim(word))
                !----------------------------------
                  case('HYDRAULIC_APERTURE_VALUE')
                    call InputReadDouble(input,option,new_fracfam%hap0)
                    call InputErrorMsg(input,option, &
                                       'HYDRAULIC_APERTURE_VALUE',error_string)
                !----------------------------------
                  case('HYDRAULIC_APERTURE_MAX')
                    call InputReadDouble(input,option,new_fracfam%hap0_max)
                    call InputErrorMsg(input,option, &
                                       'HYDRAULIC_APERTURE_MAX',error_string)
                !----------------------------------
                  case('HYDRAULIC_APERTURE_STDEV')
                    call InputReadDouble(input,option,new_fracfam%hap0_stdev)
                    call InputErrorMsg(input,option, &
                                       'HYDRAULIC_APERTURE_STDEV',error_string)
                !-----------------------------
                  case('HYDRAULIC_APERTURE_SEED')
                    call InputReadInt(input,option,new_fracfam%hap0_seed)
                    call InputErrorMsg(input,option,'HYDRAULIC_APERTURE_SEED', &
                                       error_string)
                !-----------------------------------
                  case default
                    call InputKeywordUnrecognized(input,word, &
                                                  error_string,option)
                !-----------------------------------
                end select
              enddo
            call InputPopBlock(input,option)
        !-----------------------------
          case('CENTER')
            call InputPushBlock(input,option)
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                call InputReadCard(input,option,word)
                call InputErrorMsg(input,option,'keyword',error_string)
                call StringToUpper(word)
                select case(trim(word))
                !----------------------------------
                  case('COORDINATE')
                    call InputReadDouble(input,option,new_fracfam%center%x)
                    call InputErrorMsg(input,option,'CENTER,X COORDINATE', &
                                       error_string)
                    call InputReadDouble(input,option,new_fracfam%center%y)
                    call InputErrorMsg(input,option,'CENTER,Y COORDINATE', &
                                       error_string)
                    call InputReadDouble(input,option,new_fracfam%center%z)
                    call InputErrorMsg(input,option,'CENTER,Z COORDINATE', &
                                       error_string)
                !----------------------------------
                  case('XCOORD_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%center_stdev%x)
                    call InputErrorMsg(input,option,'CENTER,XCOORD_STDEV', &
                                       error_string)
                !----------------------------------
                  case('YCOORD_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%center_stdev%y)
                    call InputErrorMsg(input,option,'CENTER,YCOORD_STDEV', &
                                       error_string)
                !----------------------------------
                  case('ZCOORD_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%center_stdev%z)
                    call InputErrorMsg(input,option,'CENTER,ZCOORD_STDEV', &
                                       error_string)
                !-----------------------------
                  case('CENTER_SEED')
                    call InputReadInt(input,option,new_fracfam%center_seed)
                    call InputErrorMsg(input,option,'CENTER_SEED',error_string)
                !-----------------------------------
                  case default
                    call InputKeywordUnrecognized(input,word, &
                                                  error_string,option)
                !-----------------------------------
                end select
              enddo
            call InputPopBlock(input,option)
        !-----------------------------
          case('NORMAL_VECTOR')
            call InputPushBlock(input,option)
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                call InputReadCard(input,option,word)
                call InputErrorMsg(input,option,'keyword',error_string)
                call StringToUpper(word)
                select case(trim(word))
                !----------------------------------
                  case('VECTOR_COORDINATES')
                    call InputReadDouble(input,option,new_fracfam%normal%x)
                    call InputErrorMsg(input,option, &
                      'NORMAL_VECTOR,X COORDINATE',error_string)
                    call InputReadDouble(input,option,new_fracfam%normal%y)
                    call InputErrorMsg(input,option, &
                      'NORMAL_VECTOR,Y COORDINATE',error_string)
                    call InputReadDouble(input,option,new_fracfam%normal%z)
                    call InputErrorMsg(input,option, &
                      'NORMAL_VECTOR,Z COORDINATE',error_string)
                !----------------------------------
                  case('XCOORD_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%normal_stdev%x)
                    call InputErrorMsg(input,option, &
                      'NORMAL_VECTOR,XCOORD_STDEV',error_string)
                !----------------------------------
                  case('YCOORD_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%normal_stdev%y)
                    call InputErrorMsg(input,option, &
                      'NORMAL_VECTOR,YCOORD_STDEV',error_string)
                !----------------------------------
                  case('ZCOORD_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%normal_stdev%z)
                    call InputErrorMsg(input,option, &
                      'NORMAL_VECTOR,ZCOORD_STDEV',error_string)
                !-----------------------------
                  case('NORMAL_SEED')
                    call InputReadInt(input,option,new_fracfam%normal_seed)
                    call InputErrorMsg(input,option,'NORMAL_SEED',error_string)
                !-----------------------------------
                  case default
                    call InputKeywordUnrecognized(input,word, &
                                                  error_string,option)
                !-----------------------------------
                end select
              enddo
            call InputPopBlock(input,option)
        !-----------------------------
          case('RADIUS')
            call InputPushBlock(input,option)
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                call InputReadCard(input,option,word)
                call InputErrorMsg(input,option,'keyword',error_string)
                call StringToUpper(word)
                select case(trim(word))
                !----------------------------------
                  case('RADIUS_XYZ')
                    call InputReadDouble(input,option,new_fracfam%radius%x)
                    call InputErrorMsg(input,option,'RADIUS_XYZ,X',error_string)
                    call InputReadDouble(input,option,new_fracfam%radius%y)
                    call InputErrorMsg(input,option,'RADIUS_XYZ,Y',error_string)
                    call InputReadDouble(input,option,new_fracfam%radius%z)
                    call InputErrorMsg(input,option,'RADIUS_XYZ,Z',error_string)
                !----------------------------------
                  case('RAD_X_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%radius_stdev%x)
                    call InputErrorMsg(input,option,'RADIUS_XYZ,RAD_X_STDEV', &
                                       error_string)
                !----------------------------------
                  case('RAD_Y_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%radius_stdev%y)
                    call InputErrorMsg(input,option,'RADIUS_XYZ,RAD_Y_STDEV', &
                                       error_string)
                !----------------------------------
                  case('RAD_Z_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%radius_stdev%z)
                    call InputErrorMsg(input,option,'RADIUS_XYZ,RAD_Z_STDEV', &
                                       error_string)
                !-----------------------------
                  case('RADIUS_SEED')
                    call InputReadInt(input,option,new_fracfam%radius_seed)
                    call InputErrorMsg(input,option,'RADIUS_SEED',error_string)
                !-----------------------------------
                  case default
                    call InputKeywordUnrecognized(input,word, &
                                                  error_string,option)
                !-----------------------------------
                end select
              enddo
            call InputPopBlock(input,option)
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)

      !------ error messaging ----------------------------------------
      if (Uninitialized(new_fracfam%id)) then
        option%io_buffer = 'ERROR: ID must be specified in ' // &
                           trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%intensity_fracfam)) then
        option%io_buffer = 'ERROR: FRACTURE_INTENSITY must be specified &
          &in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%hap0)) then
        option%io_buffer = 'ERROR: HYDRAULIC_APERTURE,HYDRAULIC_APERTURE_VALUE &
        &must be specified in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%hap0_max)) then
        option%io_buffer = 'ERROR: HYDRAULIC_APERTURE,HYDRAULIC_APERTURE_MAX &
        &must be specified in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%hap0_stdev)) then
        option%io_buffer = 'ERROR: HYDRAULIC_APERTURE,HYDRAULIC_APERTURE_STDEV &
        &must be specified in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%center%x) .or. &
          Uninitialized(new_fracfam%center%y) .or. &
          Uninitialized(new_fracfam%center%z)) then
        option%io_buffer = 'ERROR: CENTER coordinate must be specified in ' // &
                           trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%center_stdev%x) .or. &
          Uninitialized(new_fracfam%center_stdev%y) .or. &
          Uninitialized(new_fracfam%center_stdev%z)) then
        option%io_buffer = 'ERROR: CENTER,(X/Y/Z)COORD_STDEV must be &
          &specified in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%normal%x) .or. &
          Uninitialized(new_fracfam%normal%y) .or. &
          Uninitialized(new_fracfam%normal%z)) then
        option%io_buffer = 'ERROR: NORMAL_VECTOR coordinate must be &
                           &specified in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%normal_stdev%x) .or. &
          Uninitialized(new_fracfam%normal_stdev%y) .or. &
          Uninitialized(new_fracfam%normal_stdev%z)) then
        option%io_buffer = 'ERROR: NORMAL,(X/Y/Z)COORD_STDEV must be &
          &specified in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%radius%x) .or. &
          Uninitialized(new_fracfam%radius%y) .or. &
          Uninitialized(new_fracfam%radius%z)) then
        option%io_buffer = 'ERROR: RADIUS,RADIUS_XYZ must be &
                           &specified in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%radius_stdev%x) .or. &
          Uninitialized(new_fracfam%radius_stdev%y) .or. &
          Uninitialized(new_fracfam%radius_stdev%z)) then
        option%io_buffer = 'ERROR: RADIUS,RAD_(X/Y/Z)_STDEV must be &
          &specified in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif

      if (num_errors > 0) then
        write(option%io_buffer,*) num_errors
        option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                           &the FRACTURE_FAMILY block. See above to fix.'
        call PrintErrMsg(option)
      endif
      !------ add fracture family to list -----------------------------------
      if (.not.associated(pm_fracture%fracfam_list)) then
        pm_fracture%fracfam_list => new_fracfam
      else
        cur_fracfam => pm_fracture%fracfam_list
        do
          if (.not.associated(cur_fracfam)) exit
          if (.not.associated(cur_fracfam%next)) then
            cur_fracfam%next => new_fracfam
            added = PETSC_TRUE
          endif
          if (added) exit
          cur_fracfam => cur_fracfam%next
        enddo
      endif
      nullify(new_fracfam)
  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

end subroutine PMFracReadFracFam

! ************************************************************************** !

subroutine PMFracReadClipRegion(pm_fracture,input,option,keyword,error_string, &
                                found)
  !
  ! Reads input file parameters associated with the fracture model grid.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 05/06/2025
  !

  implicit none

  class(pm_fracture_type) :: pm_fracture
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  character(len=MAXWORDLENGTH) :: word

  error_string = trim(error_string) // ',CLIP_UNCONNECTED_FRACTURES'
  found = PETSC_TRUE

  select case(trim(keyword))
  !-------------------------------------
    case('CLIP_UNCONNECTED_FRACTURES')
      pm_fracture%clip_unconnected_fracs = PETSC_TRUE
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case(trim(word))
        !-----------------------------
          case('REGION_START')
            call InputReadWord(input,option,pm_fracture%start_region_name, &
                               PETSC_TRUE)
            call InputErrorMsg(input,option,'REGION',error_string)
                  !-----------------------------
          case('REGION_END')
            call InputReadWord(input,option,pm_fracture%end_region_name, &
                               PETSC_TRUE)
            call InputErrorMsg(input,option,'END_REGION',error_string)
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)
  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

end subroutine PMFracReadClipRegion

! ************************************************************************** !

subroutine PMFracReadPass2(input,option)
  !
  ! Reads input file parameters associated with the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXWORDLENGTH) :: card

  error_string = 'SUBSURFACE,GEOTHERMAL_FRACTURE_MODEL'

  input%ierr = INPUT_ERROR_NONE
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,card)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
      !--------------------
      case('FRACTURE')
        call InputSkipToEND(input,option,card)
      !--------------------
      case('CLIP_UNCONNECTED_FRACTURES')
        call InputSkipToEND(input,option,card)
      !--------------------
      case('FRACTURE_FAMILY')
        call InputSkipToEND(input,option,card)
        call InputSkipToEND(input,option,card)
        call InputSkipToEND(input,option,card)
        call InputSkipToEND(input,option,card)
        call InputSkipToEND(input,option,card)
      !--------------------
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine PMFracReadPass2

! ************************************************************************** !

subroutine GetLHSValues(mean,st_dev,isamples,seed,LHSvalues)
  !
  ! Gets a set of samples from a normal distribution using LHS method.
  !
  ! Author: Jenn Frederick
  ! Date: 11/26/2024
  use PFLOTRAN_Constants_module

  implicit none

  PetscReal :: mean, st_dev
  PetscInt :: isamples, seed
  PetscReal:: LHSvalues(:)

  PetscReal :: rndnum
  PetscInt :: f, k
  PetscBool :: switch
  PetscReal :: n0, n1
  PetscReal :: nu1, nu2
  PetscReal :: TWO_PI

  switch = mod(seed, 2) == 0
  TWO_PI = 2*PI

  do k=1,isamples
    if (.not.switch) then
      ! Generate two random numbers between (0,1)
      do f=1,seed
        nu1 = rnd4LHS(seed)
        nu2 = rnd4LHS(seed)
      enddo
      n0 = sqrt(-2.0*log(nu1)) * cos(TWO_PI*nu2)
      n1 = sqrt(-2.0*log(nu1)) * sin(TWO_PI*nu2)
      rndnum = n0*st_dev + mean
    else
      rndnum = n1*st_dev + mean
    endif
    LHSvalues(k) = rndnum
    switch = .not.switch
  enddo

end subroutine GetLHSValues

! ************************************************************************** !

function rnd4LHS(iseed)

  implicit none

  PetscInt :: iseed
  PetscReal :: rnd4LHS

  iseed = iseed*125
  iseed = iseed - (iseed/2796203) * 2796203
  rnd4LHS   = dble(iseed)/2796203.0
  return

end function rnd4LHS

! ************************************************************************** !

subroutine PMFracGenerateFracFam(fracfam)
  !
  ! Generates fractures based on the information about a fracture family.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/20/2024
  !
  use Utility_module

  implicit none

  type(fracfam_type) :: fracfam

  type(fracture_type), pointer :: new_fracture,cur_fracture
  PetscReal, pointer :: center_x(:),center_y(:),center_z(:)
  PetscReal, pointer :: normal_x(:),normal_y(:),normal_z(:)
  PetscReal, pointer :: radius_x(:),radius_y(:),radius_z(:)
  PetscReal, pointer :: hap0(:)
  PetscBool :: added

  allocate(center_x(fracfam%intensity_fracfam))
  allocate(center_y(fracfam%intensity_fracfam))
  allocate(center_z(fracfam%intensity_fracfam))
  center_x(:) = 0.d0; center_y(:) = 0.d0; center_z(:) = 0.d0
  allocate(normal_x(fracfam%intensity_fracfam))
  allocate(normal_y(fracfam%intensity_fracfam))
  allocate(normal_z(fracfam%intensity_fracfam))
  normal_x(:) = 0.d0; normal_y(:) = 0.d0; normal_z(:) = 0.d0
  allocate(radius_x(fracfam%intensity_fracfam))
  allocate(radius_y(fracfam%intensity_fracfam))
  allocate(radius_z(fracfam%intensity_fracfam))
  radius_x(:) = 0.d0; radius_y(:) = 0.d0; radius_z(:) = 0.d0
  allocate(hap0(fracfam%intensity_fracfam))
  hap0(:) = 0.d0

  call GetLHSValues(fracfam%center%x,fracfam%center_stdev%x, &
                    fracfam%intensity_fracfam,fracfam%center_seed,center_x)
  call GetLHSValues(fracfam%center%y,fracfam%center_stdev%y, &
                    fracfam%intensity_fracfam,fracfam%center_seed,center_y)
  call GetLHSValues(fracfam%center%z,fracfam%center_stdev%z, &
                    fracfam%intensity_fracfam,fracfam%center_seed,center_z)

  call GetLHSValues(fracfam%normal%x,fracfam%normal_stdev%x, &
                    fracfam%intensity_fracfam,fracfam%normal_seed,normal_x)
  call GetLHSValues(fracfam%normal%y,fracfam%normal_stdev%y, &
                    fracfam%intensity_fracfam,fracfam%normal_seed,normal_y)
  call GetLHSValues(fracfam%normal%z,fracfam%normal_stdev%z, &
                    fracfam%intensity_fracfam,fracfam%normal_seed,normal_z)

  call GetLHSValues(fracfam%radius%x,fracfam%radius_stdev%x, &
                    fracfam%intensity_fracfam,fracfam%radius_seed,radius_x)
  call GetLHSValues(fracfam%radius%y,fracfam%radius_stdev%y, &
                    fracfam%intensity_fracfam,fracfam%radius_seed,radius_y)
  call GetLHSValues(fracfam%radius%z,fracfam%radius_stdev%z, &
                    fracfam%intensity_fracfam,fracfam%radius_seed,radius_z)

  call GetLHSValues(fracfam%hap0,fracfam%hap0_stdev, &
                    fracfam%intensity_fracfam,fracfam%hap0_seed,hap0)

  do while (fracfam%nfrac_in_fam < fracfam%intensity_fracfam)
    added = PETSC_FALSE
    allocate(new_fracture)
    fracfam%nfrac_in_fam = fracfam%nfrac_in_fam + 1
    new_fracture => FractureCreate()
    new_fracture%id = fracfam%id*100 + (fracfam%nfrac_in_fam)

    new_fracture%center%x = center_x(fracfam%nfrac_in_fam)
    new_fracture%center%y = center_y(fracfam%nfrac_in_fam)
    new_fracture%center%z = center_z(fracfam%nfrac_in_fam)

    new_fracture%normal%x = normal_x(fracfam%nfrac_in_fam)
    new_fracture%normal%y = normal_y(fracfam%nfrac_in_fam)
    new_fracture%normal%z = normal_z(fracfam%nfrac_in_fam)

    new_fracture%radx = radius_x(fracfam%nfrac_in_fam)
    new_fracture%rady = radius_y(fracfam%nfrac_in_fam)
    new_fracture%radz = radius_z(fracfam%nfrac_in_fam)
    new_fracture%radx = abs(new_fracture%radx)  ! ensures radius is positive
    new_fracture%rady = abs(new_fracture%rady)  ! ensures radius is positive
    new_fracture%radz = abs(new_fracture%radz)  ! ensures radius is positive

    new_fracture%hap0 = hap0(fracfam%nfrac_in_fam)
    new_fracture%hap0 = abs(new_fracture%hap0)  ! ensures hap0 is positive
    new_fracture%hap0 = min(new_fracture%hap0,fracfam%hap0_max)

    !------ add fracture to list -------------------------
    if (.not.associated(fracfam%fracture_list)) then
      fracfam%fracture_list => new_fracture
    else
      cur_fracture => fracfam%fracture_list
      do
        if (.not.associated(cur_fracture)) exit
        if (.not.associated(cur_fracture%next)) then
          cur_fracture%next => new_fracture
          added = PETSC_TRUE
        endif
        if (added) exit
        cur_fracture => cur_fracture%next
      enddo
    endif

    nullify(new_fracture)
  enddo

  deallocate(center_x); deallocate(radius_x); deallocate(normal_x)
  deallocate(center_y); deallocate(radius_y); deallocate(normal_y)
  deallocate(center_z); deallocate(radius_z); deallocate(normal_z)
  deallocate(hap0)

end subroutine PMFracGenerateFracFam

! ************************************************************************** !

subroutine PMFracCalcK(L,hap,RM,kx,ky,kz)
  !
  ! Calculates cell anisotropic permeability from fracture permeability.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/15/2023
  !

  implicit none

  PetscReal :: L
  PetscReal :: hap
  PetscReal :: RM(3,3)
  PetscReal :: kx,ky,kz

  PetscReal :: K_frac(3,3)
  PetscReal :: K_frac_rot(3,3)
  PetscReal :: K_domain(3,3)

  K_frac(:,:) = 0.d0
  K_frac(1,1) = (1.d0/12.d0)*(hap**3.0d0)/L ! kxx
  K_frac(2,2) = (1.d0/12.d0)*(hap**3.0d0)/L ! kyy
  ! note: zero in kzz

  K_frac_rot = MATMUL(RM,K_frac)
  K_domain = MATMUL(K_frac_rot,TRANSPOSE(RM))

  kx = K_domain(1,1)
  ky = K_domain(2,2)
  kz = K_domain(3,3)

end subroutine PMFracCalcK

! ************************************************************************** !

subroutine PMFracSolve(this,time,ierr)
  !
  ! Main solve step for the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !
  use Material_Aux_module
  use TH_Aux_module
  use General_Aux_module
  use Field_module
  use Grid_module

  implicit none

  class(pm_fracture_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  type(material_auxvar_type), pointer :: material_auxvar
  type(th_auxvar_type), pointer :: th_auxvar
  type(general_auxvar_type), pointer :: general_auxvar
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(fracture_type), pointer :: cur_fracture
  PetscReal, pointer :: perm0_xx_p(:), perm0_yy_p(:), perm0_zz_p(:)
  PetscReal :: cur_temperature,L
  PetscReal :: kx,ky,kz
  PetscInt, pointer :: cellid_list(:)
  PetscInt :: icell,k

  option => this%option
  field => this%realization%field
  grid => this%realization%patch%grid

  call VecGetArrayRead(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)
  this%max_frac_kx = 0.d0; this%max_frac_ky = 0.d0; this%max_frac_kz = 0.d0
  this%sum_frac_kx = 0.d0; this%sum_frac_ky = 0.d0; this%sum_frac_kz = 0.d0
  this%max_frac_hap = 0.d0; this%sum_frac_hap = 0.d0

  cur_fracture => this%fracture_list
  do
    if (.not.associated(cur_fracture)) exit

    do k = 1,cur_fracture%ncells
      icell = cur_fracture%cell_ids(k)
      material_auxvar => &
        this%realization%patch%aux%material%auxvars(grid%nL2G(icell))
      if (option%iflowmode == G_MODE) then
        general_auxvar => &
          this%realization%patch%aux%General%auxvars(1,grid%nL2G(icell))
        cur_temperature = general_auxvar%temp
      elseif (option%iflowmode == TH_MODE) then
        th_auxvar => &
          this%realization%patch%aux%TH%auxvars(grid%nL2G(icell))
        cur_temperature = th_auxvar%temp
      endif

      if (cur_fracture%prev_temperature(k) == 0.d0) then
        ! for whatever reason, in general_mode the auxvar temperature is not
        ! initialized before the fracture pm calls initialize run.
        cur_fracture%dT(k) = 0.d0
      else
        cur_fracture%dT(k) = cur_temperature - cur_fracture%prev_temperature(k)
      endif

      ! if cells aren't too pancaked, L should be a good estimate of length
      L = (material_auxvar%volume)**(1.d0/3.d0)
      cur_fracture%dL(k) = L*cur_fracture%dT(k)*this%t_coeff

      ! update the aperture with dL and calculate new perm
      cur_fracture%hap(k) = cur_fracture%hap(k) - cur_fracture%dL(k)
      ! make sure aperature doesn't go negative
      cur_fracture%hap(k) = max(0.d0,cur_fracture%hap(k))

      ! get anisotropic domain permeability from rotation transformation
      call PMFracCalcK(L,cur_fracture%hap(k),cur_fracture%RM,kx,ky,kz)

      ! record the maximum permeability encountered so far
      this%max_frac_kx(icell) = max(kx,this%max_frac_kx(icell))
      this%max_frac_ky(icell) = max(ky,this%max_frac_ky(icell))
      this%max_frac_kz(icell) = max(kz,this%max_frac_kz(icell))
      this%max_frac_hap(icell) = &
        max(cur_fracture%hap(k),this%max_frac_hap(icell))

      ! record the sum of permeability encountered so far
      this%sum_frac_kx(icell) = this%sum_frac_kx(icell) + kx
      this%sum_frac_ky(icell) = this%sum_frac_ky(icell) + ky
      this%sum_frac_kz(icell) = this%sum_frac_kz(icell) + kz
      this%sum_frac_hap(icell) = this%sum_frac_hap(icell) + cur_fracture%hap(k)

      ! reset previous temperature for next time step
      cur_fracture%prev_temperature(k) = cur_temperature
    enddo

    cur_fracture => cur_fracture%next
  enddo

  if (this%clip_unconnected_fracs) then
    cellid_list => this%unique_connected_cell_ids ! local cell id list
  else
    cellid_list => this%allfrac_unique_cell_ids ! local cell id list
  endif

  ! update domain material permeability of cells with fractures
  if (associated(cellid_list)) then
    do k = 1,size(cellid_list)
      icell = cellid_list(k)
      material_auxvar => &
        this%realization%patch%aux%material%auxvars(grid%nL2G(icell))
      if (this%frac_intersection_type == 2) then ! maximum
        perm0_xx_p(icell) = perm0_xx_p(icell) + this%max_frac_kx(icell)
        perm0_yy_p(icell) = perm0_yy_p(icell) + this%max_frac_ky(icell)
        perm0_zz_p(icell) = perm0_zz_p(icell) + this%max_frac_kz(icell)
        material_auxvar%fracture_aperture = this%max_frac_hap(icell)
      endif
      if (this%frac_intersection_type == 1) then ! summation
        perm0_xx_p(icell) = perm0_xx_p(icell) + this%sum_frac_kx(icell)
        perm0_yy_p(icell) = perm0_yy_p(icell) + this%sum_frac_ky(icell)
        perm0_zz_p(icell) = perm0_zz_p(icell) + this%sum_frac_kz(icell)
        material_auxvar%fracture_aperture = this%sum_frac_hap(icell)
      endif
    enddo
  endif

  if (this%update_material_auxvar_perm) then
    if (associated(cellid_list)) then
      do k = 1,size(cellid_list)
        icell = cellid_list(k)
        material_auxvar => &
          this%realization%patch%aux%material%auxvars(grid%nL2G(icell))
        material_auxvar%permeability(1) = perm0_xx_p(icell)
        material_auxvar%permeability(2) = perm0_yy_p(icell)
        material_auxvar%permeability(3) = perm0_zz_p(icell)
      enddo
    endif
  endif

  call VecRestoreArrayRead(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)

  call PMFracPrintStats(this)
  
  ierr = 0 ! If this is not set to zero, TS_STOP_FAILURE occurs

end subroutine PMFracSolve

! ************************************************************************** !

subroutine PMFracPrintStats(this)
  !
  ! Adds fracture network statistics to .info file.
  !
  ! Author: Jenn Frederick
  ! Date: 02/04/2026
  !
  use Material_Aux_module
  use Grid_module
  use Utility_module

  class(pm_fracture_type) :: this

  character(len=MAXWORDLENGTH) :: filename
  type(material_auxvar_type), pointer :: material_auxvar
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(grid_type), pointer :: grid
  PetscReal, pointer :: frac_aperture_lcl(:),frac_aperture_glb(:)
  PetscReal, pointer :: frac_perm_x_lcl(:),frac_perm_x_glb(:)
  PetscReal, pointer :: frac_perm_y_lcl(:),frac_perm_y_glb(:)
  PetscReal, pointer :: frac_perm_z_lcl(:),frac_perm_z_glb(:)
  PetscReal, pointer :: fracture_volume_glb(:)
  PetscReal :: fracture_volume_lcl,L2
  PetscReal :: min_frac_aperture,max_frac_aperture,mean_frac_aperture
  PetscReal :: range_frac_aperture,midrange_frac_aperture
  PetscReal :: variance_frac_aperture
  PetscReal :: skewness_frac_aperture_p, skewness_frac_aperture_a
  PetscReal :: perm_sum_x,perm_sum_y,perm_sum_z
  PetscReal :: min_perm_x,max_perm_x
  PetscReal :: min_perm_y,max_perm_y
  PetscReal :: min_perm_z,max_perm_z
  PetscReal :: range_perm_x,range_perm_y,range_perm_z
  PetscReal :: midrange_perm_x,midrange_perm_y,midrange_perm_z
  PetscReal :: mean_perm_x,mean_perm_y,mean_perm_z
  PetscReal :: variance_perm_x,variance_perm_y,variance_perm_z
  PetscReal :: skewness_perm_x_p, skewness_perm_x_a
  PetscReal :: skewness_perm_y_p, skewness_perm_y_a
  PetscReal :: skewness_perm_z_p, skewness_perm_z_a
  PetscInt, pointer :: lcl_cells(:),cell_positions(:)
  PetscInt :: i,k,fid,nlcl_cells,nglb_cells
  PetscErrorCode :: ierr

  if (this%clip_unconnected_fracs) then
    nlcl_cells = size(this%unique_connected_cell_ids) ! local cell ids
    nglb_cells = size(this%unique_connected_cell_ids_glb) ! global cell ids
  else
    nlcl_cells = size(this%allfrac_unique_cell_ids)
    nglb_cells = size(this%allfrac_unique_cell_ids_glb)
  endif
  allocate(frac_aperture_lcl(nlcl_cells)) ! ok if nlcl_cells = 0
  allocate(frac_aperture_glb(nglb_cells)) ! ok if nlcl_cells = 0
  allocate(frac_perm_x_lcl(nlcl_cells))
  allocate(frac_perm_y_lcl(nlcl_cells))
  allocate(frac_perm_z_lcl(nlcl_cells))
  allocate(frac_perm_x_glb(nglb_cells))
  allocate(frac_perm_y_glb(nglb_cells))
  allocate(frac_perm_z_glb(nglb_cells))

  ! filename needs to be consistant w/ filename located at end of PMFracSetup()
  filename = trim(this%option%global_prefix)//'.info'
  fid = 19181723
110 format(100es16.8)
  OPEN(UNIT=fid,FILE=filename,STATUS='OLD',POSITION='APPEND',ACTION='WRITE')

  fracture_volume_lcl = 0.d0
  min_perm_x = 1.d20; min_perm_y = 1.d20; min_perm_z = 1.d20
  max_perm_x = 0.d0; max_perm_y = 0.d0; max_perm_z = 0.d0
  perm_sum_x = 0.d0; perm_sum_y = 0.d0; perm_sum_z = 0.d0

  material_auxvars => this%realization%patch%aux%material%auxvars
  grid => this%realization%patch%grid

  do k = 1,nlcl_cells
    if (this%clip_unconnected_fracs) then
      material_auxvar => &
               material_auxvars(grid%nL2G(this%unique_connected_cell_ids(k)))
    else
      material_auxvar => &
               material_auxvars(grid%nL2G(this%allfrac_unique_cell_ids(k)))
    endif
    L2 = ((material_auxvar%volume)**(1.d0/3.d0))**(2)! [m^2]
    frac_aperture_lcl(k) = material_auxvar%fracture_aperture
    frac_perm_x_lcl(k) = material_auxvar%permeability(1)
    frac_perm_y_lcl(k) = material_auxvar%permeability(2)
    frac_perm_z_lcl(k) = material_auxvar%permeability(3)
    fracture_volume_lcl = fracture_volume_lcl + (L2*frac_aperture_lcl(k)) ![m^3]
  enddo

  ! create global arrays of aperture and perm local vectors
  allocate(fracture_volume_glb(this%option%comm%size))
  if (nglb_cells > 0) then
    allocate(lcl_cells(this%option%comm%size))
    call MPI_Allgather(nlcl_cells,ONE_INTEGER,MPI_INTEGER,lcl_cells, &
      ONE_INTEGER,MPI_INTEGER,this%option%mycomm,ierr);CHKERRQ(ierr)
    call MPI_Barrier(this%option%mycomm,ierr);CHKERRQ(ierr)
    allocate(cell_positions(this%option%comm%size))
    do i = 1,this%option%comm%size
      if (i == 1) then
        cell_positions(i) = 0
      else
        cell_positions(i) = cell_positions(i-1) + lcl_cells(i-1)
      endif
    enddo
    ! fracture aperture -------------------------------------------------------
    call MPI_Gatherv(frac_aperture_lcl,nlcl_cells,MPI_DOUBLE_PRECISION, &
      frac_aperture_glb,lcl_cells,cell_positions,MPI_DOUBLE_PRECISION, &
      ZERO_INTEGER,this%option%mycomm,ierr);CHKERRQ(ierr)
         ! only rank zero knows frac_aperture_glb, so we need to share it
    call MPI_Barrier(this%option%mycomm,ierr);CHKERRQ(ierr)
    call MPI_Bcast(frac_aperture_glb,nglb_cells,MPI_DOUBLE_PRECISION, &
                   ZERO_INTEGER,this%option%mycomm,ierr);CHKERRQ(ierr)
    ! fracture x-permeability -------------------------------------------------
    call MPI_Gatherv(frac_perm_x_lcl,nlcl_cells,MPI_DOUBLE_PRECISION, &
      frac_perm_x_glb,lcl_cells,cell_positions,MPI_DOUBLE_PRECISION, &
      ZERO_INTEGER,this%option%mycomm,ierr);CHKERRQ(ierr)
         ! only rank zero knows frac_perm_x_glb, so we need to share it
    call MPI_Barrier(this%option%mycomm,ierr);CHKERRQ(ierr)
    call MPI_Bcast(frac_perm_x_glb,nglb_cells,MPI_DOUBLE_PRECISION, &
                   ZERO_INTEGER,this%option%mycomm,ierr);CHKERRQ(ierr)
    ! fracture y-permeability -------------------------------------------------
    call MPI_Gatherv(frac_perm_y_lcl,nlcl_cells,MPI_DOUBLE_PRECISION, &
      frac_perm_y_glb,lcl_cells,cell_positions,MPI_DOUBLE_PRECISION, &
      ZERO_INTEGER,this%option%mycomm,ierr);CHKERRQ(ierr)
         ! only rank zero knows frac_perm_y_glb, so we need to share it
    call MPI_Barrier(this%option%mycomm,ierr);CHKERRQ(ierr)
    call MPI_Bcast(frac_perm_y_glb,nglb_cells,MPI_DOUBLE_PRECISION, &
                   ZERO_INTEGER,this%option%mycomm,ierr);CHKERRQ(ierr)
    ! fracture z-permeability -------------------------------------------------
    call MPI_Gatherv(frac_perm_z_lcl,nlcl_cells,MPI_DOUBLE_PRECISION, &
      frac_perm_z_glb,lcl_cells,cell_positions,MPI_DOUBLE_PRECISION, &
      ZERO_INTEGER,this%option%mycomm,ierr);CHKERRQ(ierr)
         ! only rank zero knows frac_perm_z_glb, so we need to share it
    call MPI_Barrier(this%option%mycomm,ierr);CHKERRQ(ierr)
    call MPI_Bcast(frac_perm_z_glb,nglb_cells,MPI_DOUBLE_PRECISION, &
                   ZERO_INTEGER,this%option%mycomm,ierr);CHKERRQ(ierr)
    ! fracture volume ---------------------------------------------------------
    call MPI_Allgather(fracture_volume_lcl,ONE_INTEGER,MPI_DOUBLE_PRECISION, &
      fracture_volume_glb,ONE_INTEGER,MPI_DOUBLE_PRECISION, &
      this%option%mycomm,ierr);CHKERRQ(ierr)
    call MPI_Barrier(this%option%mycomm,ierr);CHKERRQ(ierr)
    !call MPI_Barrier(this%option%mycomm,ierr);CHKERRQ(ierr)
    !call MPI_Allreduce(MPI_IN_PLACE,fracture_volume,this%option%comm%size, &
    !  MPI_DOUBLE_PRECISION,MPI_SUM,this%option%mycomm,ierr);CHKERRQ(ierr)
  endif

  perm_sum_x = SUM(frac_perm_x_glb)
  perm_sum_y = SUM(frac_perm_y_glb)
  perm_sum_z = SUM(frac_perm_z_glb)

  max_frac_aperture = maxval(frac_aperture_glb)
  min_frac_aperture = minval(frac_aperture_glb)
  max_perm_x = maxval(frac_perm_x_glb); min_perm_x = minval(frac_perm_x_glb)
  max_perm_y = maxval(frac_perm_y_glb); min_perm_y = minval(frac_perm_y_glb)
  max_perm_z = maxval(frac_perm_z_glb); min_perm_z = minval(frac_perm_z_glb)

  range_frac_aperture = max_frac_aperture-min_frac_aperture
  midrange_frac_aperture = range_frac_aperture/2.d0

  range_perm_x = max_perm_x-min_perm_x
  midrange_perm_x = range_perm_x/2.d0
  range_perm_y = max_perm_y-min_perm_y
  midrange_perm_y = range_perm_y/2.d0
  range_perm_z = max_perm_z-min_perm_z
  midrange_perm_z = range_perm_z/2.d0

  mean_frac_aperture = SUM(frac_aperture_glb)/nglb_cells
  mean_perm_x = perm_sum_x/nglb_cells
  mean_perm_y = perm_sum_y/nglb_cells
  mean_perm_z = perm_sum_z/nglb_cells

  variance_frac_aperture = &
                  SUM((frac_aperture_glb-mean_frac_aperture)**2)/(nglb_cells-1)
  variance_perm_x = SUM((frac_perm_x_glb-mean_perm_x)**2)/(nglb_cells-1)
  variance_perm_y = SUM((frac_perm_y_glb-mean_perm_y)**2)/(nglb_cells-1)
  variance_perm_z = SUM((frac_perm_z_glb-mean_perm_z)**2)/(nglb_cells-1)

  skewness_frac_aperture_p = Skewness(frac_aperture_glb)
  skewness_frac_aperture_a = Skewness(frac_aperture_glb,PETSC_TRUE)
  skewness_perm_x_p = Skewness(frac_perm_x_glb)
  skewness_perm_x_a = Skewness(frac_perm_x_glb,PETSC_TRUE)
  skewness_perm_y_p = Skewness(frac_perm_y_glb)
  skewness_perm_y_a = Skewness(frac_perm_y_glb,PETSC_TRUE)
  skewness_perm_z_p = Skewness(frac_perm_z_glb)
  skewness_perm_z_a = Skewness(frac_perm_z_glb,PETSC_TRUE)

  ! header for these values is located at end of PMFracSetup()
  ! TIME [tunit];
  write(fid,110,advance="no") this%option%time/this%output_option%tconv
  write(fid,'(a)',advance="no") ","
  ! VOL [m^3] OF FRACTURES;
  write(fid,110,advance="no") SUM(fracture_volume_glb)
  write(fid,'(a)',advance="no") ","
  ! MIN FRAC APERTURE [m];
  write(fid,110,advance="no") min_frac_aperture
  write(fid,'(a)',advance="no") ","
  ! MAX FRAC APERTURE [m];
  write(fid,110,advance="no") max_frac_aperture
  write(fid,'(a)',advance="no") ","
  ! RANGE FRAC APERTURE [m];
  write(fid,110,advance="no") range_frac_aperture
  write(fid,'(a)',advance="no") ","
  ! MIDRANGE FRAC APERTURE [m];
  write(fid,110,advance="no") midrange_frac_aperture
  write(fid,'(a)',advance="no") ","
  ! MEAN FRAC APERTURE [m];
  write(fid,110,advance="no") mean_frac_aperture
  write(fid,'(a)',advance="no") ","
  ! VARIANCE FRAC APERTURE [m];
  write(fid,110,advance="no") variance_frac_aperture
  write(fid,'(a)',advance="no") ","
  ! POP SKEWNESS FRAC APERTURE [m];
  write(fid,110,advance="no") skewness_frac_aperture_p
  write(fid,'(a)',advance="no") ","
  ! SAMPLE SKEWNESS FRAC APERTURE [m];
  write(fid,110,advance="no") skewness_frac_aperture_a
  write(fid,'(a)',advance="no") ","
  ! MIN FRAC PERM-X [m];
  write(fid,110,advance="no") min_perm_x
  write(fid,'(a)',advance="no") ","
  ! MAX FRAC PERM-X [m];
  write(fid,110,advance="no") max_perm_x
  write(fid,'(a)',advance="no") ","
  ! RANGE FRAC PERM-X [m];
  write(fid,110,advance="no") range_perm_x
  write(fid,'(a)',advance="no") ","
  ! MIDRANGE FRAC PERM-X [m];
  write(fid,110,advance="no") midrange_perm_x
  write(fid,'(a)',advance="no") ","
  ! MEAN FRAC PERM-X [m];
  write(fid,110,advance="no") mean_perm_x
  write(fid,'(a)',advance="no") ","
  ! VARIANCE FRAC PERM-X [m];
  write(fid,110,advance="no") variance_perm_x
  write(fid,'(a)',advance="no") ","
  ! POP SKEWNESS FRAC PERM-X [m];
  write(fid,110,advance="no") skewness_perm_x_p
  write(fid,'(a)',advance="no") ","
  ! SAMPLE SKEWNESS FRAC PERM-X [m];
  write(fid,110,advance="no") skewness_perm_x_a
  write(fid,'(a)',advance="no") ","
  ! MIN FRAC PERM-Y [m];
  write(fid,110,advance="no") min_perm_y
  write(fid,'(a)',advance="no") ","
  ! MAX FRAC PERM-Y [m];
  write(fid,110,advance="no") max_perm_y
  write(fid,'(a)',advance="no") ","
  ! RANGE FRAC PERM-Y [m];
  write(fid,110,advance="no") range_perm_y
  write(fid,'(a)',advance="no") ","
  ! MIDRANGE FRAC PERM-Y [m];
  write(fid,110,advance="no") midrange_perm_y
  write(fid,'(a)',advance="no") ","
  ! MEAN FRAC PERM-Y [m];
  write(fid,110,advance="no") mean_perm_y
  write(fid,'(a)',advance="no") ","
  ! VARIANCE FRAC PERM-Y [m];
  write(fid,110,advance="no") variance_perm_y
  write(fid,'(a)',advance="no") ","
  ! POP SKEWNESS FRAC PERM-Y [m];
  write(fid,110,advance="no") skewness_perm_y_p
  write(fid,'(a)',advance="no") ","
  ! SAMPLE SKEWNESS FRAC PERM-Y [m];
  write(fid,110,advance="no") skewness_perm_y_a
  write(fid,'(a)',advance="no") ","
  ! MIN FRAC PERM-Z [m];
  write(fid,110,advance="no") min_perm_z
  write(fid,'(a)',advance="no") ","
  ! MAX FRAC PERM-Z [m];
  write(fid,110,advance="no") max_perm_z
  write(fid,'(a)',advance="no") ","
  ! RANGE FRAC PERM-Z [m];
  write(fid,110,advance="no") range_perm_z
  write(fid,'(a)',advance="no") ","
  ! MIDRANGE FRAC PERM-Z [m];
  write(fid,110,advance="no") midrange_perm_z
  write(fid,'(a)',advance="no") ","
  ! MEAN FRAC PERM-Z [m];
  write(fid,110,advance="no") mean_perm_z
  write(fid,'(a)',advance="no") ","
  ! VARIANCE FRAC PERM-Z [m];
  write(fid,110,advance="no") variance_perm_z
  write(fid,'(a)',advance="no") ","
  ! POP SKEWNESS FRAC PERM-Z [m];
  write(fid,110,advance="no") skewness_perm_z_p
  write(fid,'(a)',advance="no") ","
  ! SAMPLE SKEWNESS FRAC PERM-Z [m];
  write(fid,110,advance="yes") skewness_perm_z_a
  ! make sure last write() says advance="yes", all others "no"

  close(fid)

  deallocate(fracture_volume_glb)
  deallocate(lcl_cells); deallocate(cell_positions)
  deallocate(frac_perm_z_glb); deallocate(frac_perm_y_glb)
  deallocate(frac_perm_x_glb)
  deallocate(frac_perm_z_lcl); deallocate(frac_perm_y_lcl)
  deallocate(frac_perm_x_lcl)
  deallocate(frac_aperture_glb)
  deallocate(frac_aperture_lcl)

end subroutine PMFracPrintStats
! ************************************************************************** !
subroutine PMFracGeomPrintStats(this)
  !
  ! Adds fracture netowrk geometry information to .frac_geom_info file.
  !
  ! Author: Elyce Bayat
  ! Date: 02/25/2026
  !
  use Material_Aux_module
  use Grid_module
  use Utility_module

  class(pm_fracture_type) :: this

  character(len=MAXWORDLENGTH) :: filename_frac_info ! EJB
  type(fracture_type), pointer :: cur_fracture
  PetscInt :: fid_frac_info, frac_count

  ! filename_frac_info needs to be consistent w/ filename_frac_info located at end of PMFracSetup()
  ! EJB: filename for printing out fracture geometry information: 
  filename_frac_info = trim(this%option%global_prefix)//'.frac_geom_info'
  fid_frac_info = 20846713 ! fid here also needs to be consistent w one at end of PMFracSetup()
110 format(100es16.8)  
  OPEN(UNIT=fid_frac_info,FILE=filename_frac_info,STATUS='OLD',POSITION='APPEND',ACTION='WRITE')

  ! Loop over each fracture in the domain and print its properties
  frac_count = 1
  cur_fracture => this%fracture_list
  do 
    if (.not.associated(cur_fracture)) exit

    ! data to write is: 
    ! fracture_id, cx, cy, cz, nx, ny, nz, rx, ry, rz, aperture
    write(fid_frac_info,'(I8)',advance="no") cur_fracture%id
    write(fid_frac_info,'(a)', advance="no") ","
    
    ! center_x, center_y, center_z
    write(fid_frac_info,110,advance="no") cur_fracture%center%x
    write(fid_frac_info,'(a)', advance="no") ","
    write(fid_frac_info,110,advance="no") cur_fracture%center%y
    write(fid_frac_info,'(a)', advance="no") ","
    write(fid_frac_info,110,advance="no") cur_fracture%center%z
    write(fid_frac_info,'(a)', advance="no") ","

    ! normal_x, normal_y, normal_z
    write(fid_frac_info,110,advance="no") cur_fracture%normal%x
    write(fid_frac_info,'(a)', advance="no") ","
    write(fid_frac_info,110,advance="no") cur_fracture%normal%y
    write(fid_frac_info,'(a)', advance="no") ","
    write(fid_frac_info,110,advance="no") cur_fracture%normal%z
    write(fid_frac_info,'(a)', advance="no") ","

    ! radius_x, radius_y, radius_z
    write(fid_frac_info,110,advance="no") cur_fracture%radx
    write(fid_frac_info,'(a)', advance="no") ","
    write(fid_frac_info,110,advance="no") cur_fracture%rady
    write(fid_frac_info,'(a)', advance="no") ","
    write(fid_frac_info,110,advance="no") cur_fracture%radz
    write(fid_frac_info,'(a)', advance="no") ","

    ! initial aperture    
    write(fid_frac_info,110,advance="no") cur_fracture%hap0
    write(fid_frac_info,'(a)', advance="no") ","

    ! strike
    write(fid_frac_info, 110, advance="no") cur_fracture%strike%x
    write(fid_frac_info,'(a)', advance="no") ","
    write(fid_frac_info, 110, advance="no") cur_fracture%strike%y
    write(fid_frac_info,'(a)', advance="no") ","
    write(fid_frac_info, 110, advance="no") cur_fracture%strike%z
    write(fid_frac_info,'(a)', advance="no") ","

    ! dip
    write(fid_frac_info, 110, advance="yes") cur_fracture%dip
    
    this%option%io_buffer = 'Frac count: ' // StringWrite(frac_count)
    call PrintMsg(this%option)

    frac_count = frac_count + 1
    cur_fracture => cur_fracture%next

  enddo

  CLOSE(fid_frac_info)
end subroutine PMFracGeomPrintStats

! ************************************************************************** !

subroutine PMFracSetPlotVariables(list,this)
  !
  ! Adds variables to be printed for plotting.
  !
  ! Author: Jenn Frederick
  ! Date: 04/23/2025
  !
  use Output_Aux_module
  use Variables_module

  type(output_variable_list_type), pointer :: list
  class(pm_fracture_type) :: this

  character(len=MAXWORDLENGTH) :: name, units

  name = 'Fracture Aperture'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units,FRACTURE_APERTURE)

end subroutine PMFracSetPlotVariables

! ************************************************************************** !

subroutine PMFracDestroy(this)
  !
  ! Destroys objects in the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/18/2023
  !
  use Utility_module, only : DeallocateArray

  implicit none

  class(pm_fracture_type) :: this

  type(fracture_type), pointer :: cur_fracture

  if (associated(this%allfrac_cell_ids)) then
    call DeallocateArray(this%allfrac_cell_ids)
  endif
  if (associated(this%allfrac_unique_cell_ids)) then
    call DeallocateArray(this%allfrac_unique_cell_ids)
  endif
  if (associated(this%allfrac_unique_cell_ids_glb)) then
    call DeallocateArray(this%allfrac_unique_cell_ids_glb)
  endif
  if (associated(this%frac_common_cell_ids)) then
    call DeallocateArray(this%frac_common_cell_ids)
  endif
  if (associated(this%connected_to_start)) then
    call DeallocateArray(this%connected_to_start)
  endif
  if (associated(this%connected_to_end)) then
    call DeallocateArray(this%connected_to_end)
  endif
  if (associated(this%unique_connected_cell_ids_glb)) then
    call DeallocateArray(this%unique_connected_cell_ids_glb)
  endif
  if (associated(this%clip_region_cell_ids)) then
    call DeallocateArray(this%clip_region_cell_ids)
  endif
  if (associated(this%max_frac_kx)) then
    call DeallocateArray(this%max_frac_kx)
  endif
  if (associated(this%max_frac_ky)) then
    call DeallocateArray(this%max_frac_ky)
  endif
  if (associated(this%max_frac_ky)) then
    call DeallocateArray(this%max_frac_kz)
  endif
  if (associated(this%sum_frac_kx)) then
    call DeallocateArray(this%sum_frac_kx)
  endif
  if (associated(this%sum_frac_ky)) then
    call DeallocateArray(this%sum_frac_ky)
  endif
  if (associated(this%sum_frac_ky)) then
    call DeallocateArray(this%sum_frac_kz)
  endif

  cur_fracture => this%fracture_list
  do
    if (.not.associated(cur_fracture)) exit

    if (cur_fracture%ncells > 0) call PMFracStripFrac(cur_fracture)

    cur_fracture => cur_fracture%next
  enddo

end subroutine PMFracDestroy

! ************************************************************************** !

subroutine PMFracStripFrac(frac)
  !
  ! Destroys objects in the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/18/2023
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(fracture_type), pointer :: frac

  call DeallocateArray(frac%kx0)
  call DeallocateArray(frac%ky0)
  call DeallocateArray(frac%kx0)
  call DeallocateArray(frac%hap)
  call DeallocateArray(frac%cell_ids)
  call DeallocateArray(frac%prev_temperature)
  call DeallocateArray(frac%dT)
  call DeallocateArray(frac%dL)

end subroutine PMFracStripFrac

! ************************************************************************** !

end module PM_Fracture_class
