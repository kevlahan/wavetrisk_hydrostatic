program mpi_identity
  use mpi
  implicit none
  character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: version
  integer :: ierr,n
  call MPI_Init(ierr)
  if(ierr/=MPI_SUCCESS) error stop 'MPI identity initialization failed'
  call MPI_Get_library_version(version,n,ierr)
  if(ierr/=MPI_SUCCESS.or.n<1) error stop 'MPI library identity unavailable'
  write(*,'(a)') version(:n)
  call MPI_Finalize(ierr)
  if(ierr/=MPI_SUCCESS) error stop 'MPI identity finalization failed'
end program
