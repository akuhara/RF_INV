program main
  use params
  use mt19937
  use fftw
  use model
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
  
  ! Set verbose mode for rank 0
  verb = .false.
  if (rank == 0) verb = .true.  

  ! Read parameters from file
  call get_params(verb, "params.in") ! if rank == 0 -> verbose
  
  ! Read observed files
  call read_obs(verb)
  
  ! Initialize random number generator
  call sgrnd(iseed)

  ! Initialize FFTW
  call init_fftw()
  
  ! Read reference velocity model
  call read_ref_model(verb)

  ! Finish
  call mpi_finalize(ierr)
  
  stop
end program main
