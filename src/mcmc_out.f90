module mcmc_out
  use params
  use pt_mcmc
  implicit none 
  
contains
  
  subroutine output_results(rank, verb)
    implicit none 
    include "mpif.h"
    integer, intent(in) :: rank
    logical, intent(in) :: verb
    integer :: nk_sum(k_max), ierr, ik
    
    call mpi_reduce(nk, nk_sum, k_max, MPI_INTEGER4, MPI_SUM, &
         & 0, MPI_COMM_WORLD, ierr)

    if (rank == 0) then
       open(io_nk, file = 'num_interface.ppd', status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0,*)"ERROR: cannot create num_interface.ppd"
          call mpi_finalize(ierr)
          stop
       end if
       do ik = 1, k_max - 1
          write(io_nk,*) ik, nk_sum(ik)
       end do
       close(io_nk)
    end if

    return 
  end subroutine output_results
  
end module mcmc_out
