program main
  use params
  implicit none 
  include "mpif.h"
  integer :: nproc, rank, ierr
  logical :: verb

  ! Initialize MPI 
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank,  ierr)
  if (nproc < 2) then
     write(0,*)"ERROR: at least 2 processor is necesarry!"
     call mpi_finalize(ierr)
     stop
  end if
  
  ! Read parameters from file
  verb = .false.
  if (rank == 0) verb = .true.  
  call get_params(verb, "params.in") ! if rank == 0 -> verbose
  
  
  

  
  ! Finish
  call mpi_finalize(ierr)
  
  stop
end program main
