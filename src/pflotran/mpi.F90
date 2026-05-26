module mpi
  implicit none

  ! Communicators (match PETSc mpiuni convention)
  integer, parameter :: MPI_COMM_NULL  = 0
  integer, parameter :: MPI_COMM_SELF  = 1
  integer, parameter :: MPI_COMM_WORLD = 2

  ! Groups / handles
  integer, parameter :: MPI_GROUP_NULL = 0

  ! Datatypes
  integer, parameter :: MPI_INTEGER           = 0
  integer, parameter :: MPI_DOUBLE_PRECISION  = 0
  integer, parameter :: MPI_2DOUBLE_PRECISION = 0
  integer, parameter :: MPI_LOGICAL           = 0
  integer, parameter :: MPI_C_BOOL            = 0
  integer, parameter :: MPI_CHARACTER         = 0

  ! Reduction ops
  integer, parameter :: MPI_SUM    = 0
  integer, parameter :: MPI_MIN    = 0
  integer, parameter :: MPI_MAX    = 0
  integer, parameter :: MPI_LOR    = 0
  integer, parameter :: MPI_LAND   = 0
  integer, parameter :: MPI_MAXLOC = 0

  ! Wildcards / special ranks (match PETSc mpiuni)
  integer, parameter :: MPI_PROC_NULL  = -1
  integer, parameter :: MPI_ANY_SOURCE = -2
  integer, parameter :: MPI_ANY_TAG    = -1
  integer, parameter :: MPI_UNDEFINED  = -32766

  ! In-place marker
  integer, parameter :: MPI_IN_PLACE = -1

  ! Status handling (MPI_Status in PETSc mpiuni is 3 ints: source, tag, error)
  integer, parameter :: MPI_STATUS_SIZE = 3
  integer, parameter :: MPI_SOURCE = 1
  integer, parameter :: MPI_TAG    = 2
  integer, parameter :: MPI_ERROR  = 3
  integer, parameter :: MPI_STATUS_IGNORE(MPI_STATUS_SIZE) = 0

  ! Info objects
  integer, parameter :: MPI_INFO_NULL = 0

end module mpi