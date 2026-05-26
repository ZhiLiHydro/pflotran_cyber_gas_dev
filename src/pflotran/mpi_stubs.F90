subroutine MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm, ierr)
  implicit none
  ! Deliberately loose dummy declarations to avoid interface/rank/type checking.
  ! We do not touch buffers (single-rank, no-op reduction).
  integer :: sendbuf(*), recvbuf(*)
  integer :: count, datatype, op, comm, ierr
  ierr = 0
end subroutine MPI_Allreduce

subroutine MPI_Bcast(buf, count, datatype, root, comm, ierr)
  implicit none
  integer :: buf(*)
  integer :: count, datatype, root, comm, ierr
  ierr = 0
end subroutine MPI_Bcast

subroutine MPI_Comm_rank(comm, rank, ierr)
  implicit none
  integer :: comm, rank, ierr
  rank = 0
  ierr = 0
end subroutine MPI_Comm_rank

subroutine MPI_Comm_size(comm, size, ierr)
  implicit none
  integer :: comm, size, ierr
  size = 1
  ierr = 0
end subroutine MPI_Comm_size

subroutine MPI_Barrier(comm, ierr)
  implicit none
  integer :: comm, ierr
  ierr = 0
end subroutine MPI_Barrier

subroutine MPI_Init(ierr)
  implicit none
  integer :: ierr
  ierr = 0
end subroutine MPI_Init

subroutine MPI_Finalize(ierr)
  implicit none
  integer :: ierr
  ierr = 0
end subroutine MPI_Finalize

subroutine MPI_Group_incl(group, n, ranks, newgroup, ierr)
  implicit none
  integer :: group, n, newgroup, ierr
  integer :: ranks(*)
  ! Single-rank / stub behavior: return a non-error and a "valid" group handle
  newgroup = group
  ierr = 0
end subroutine MPI_Group_incl

subroutine MPI_Comm_group(comm, group, ierr)
  implicit none
  integer :: comm, group, ierr
  ! Return a non-null group for any non-null communicator
  if (comm .eq. 0) then
    group = 0
  else
    group = 1
  end if
  ierr = 0
end subroutine MPI_Comm_group

subroutine MPI_Group_free(group, ierr)
  implicit none
  integer :: group, ierr
  group = 0
  ierr = 0
end subroutine MPI_Group_free

subroutine MPI_Group_rank(group, rank, ierr)
  implicit none
  integer :: group, rank, ierr
  rank = 0
  ierr = 0
end subroutine MPI_Group_rank

subroutine MPI_Group_size(group, size, ierr)
  implicit none
  integer :: group, size, ierr
  size = 1
  ierr = 0
end subroutine MPI_Group_size

subroutine MPI_Send(buf, count, datatype, dest, tag, comm, ierr)
  implicit none
  integer :: buf(*), count, datatype, dest, tag, comm, ierr
  ierr = 0
end subroutine MPI_Send

subroutine MPI_Recv(buf, count, datatype, source, tag, comm, status, ierr)
  implicit none
  integer :: buf(*), count, datatype, source, tag, comm, ierr
  integer :: status(*)
  ierr = 0
end subroutine MPI_Recv

subroutine MPI_Abort(comm, errorcode, ierr)
  implicit none
  integer :: comm, errorcode, ierr
  ierr = 0
  stop 'MPI_Abort called'
end subroutine MPI_Abort