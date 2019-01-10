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
    integer :: nk_sum(k_max), ierr, ik, itrc, it, i
    integer :: namp_sum(nbin_amp, nsmp, ntrc), nmod_sum
    integer :: nprop_sum(ntype), naccept_sum(ntype)
    
    call mpi_reduce(nmod, nmod_sum, 1, MPI_INTEGER4, MPI_SUM, &
         & 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(nprop, nprop_sum, ntype, MPI_INTEGER4, MPI_SUM, &
         & 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(naccept, naccept_sum, ntype, MPI_INTEGER4, MPI_SUM, &
         & 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(nk, nk_sum, k_max, MPI_INTEGER4, MPI_SUM, &
         & 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(namp, namp_sum, nbin_amp * nsmp * ntrc, &
         & MPI_INTEGER4, MPI_SUM, 0, MPI_COMM_WORLD, ierr)


    if (rank == 0) then
       
       if (verb) then
          write(*,*)"--- Summary ---"
          write(*,*)"# of models to be sampled:", nmod_sum
          write(*,*)"# of birth proposal:", &
               & naccept_sum(1), "/", nprop_sum(1)
          write(*,*)"# of death proposal:", &
               & naccept_sum(2), "/", nprop_sum(2)
          write(*,*)"# of moving interface proposal:", &
               & naccept_sum(3), "/", nprop_sum(3)
          write(*,*)"# of perturbing Vp proposal:", &
               & naccept_sum(4), "/", nprop_sum(4)
          write(*,*)"# of perturbing Vs proposal:", &
               & naccept_sum(5), "/", nprop_sum(5)
          write(*,*)"# of perturbing noise sigma proposal:", &
               & naccept_sum(6), "/", nprop_sum(6)
          write(*,*)
       end if
       
       ! # of layer interfaces
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


       ! Synthetic RFs
       open(io_syn, file = 'syn_trace.ppd', status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0,*)"ERROR: cannot create syn_trace.ppd"
          call mpi_finalize(ierr)
          stop
       end if
       do itrc = 1, ntrc
          do it = 1, nsmp
             do i = 1, nbin_amp
                write(io_syn,*) (it-1) * delta + t_start, &
                     & amp_min + (i - 0.5d0) * dbin_amp, &
                     & namp_sum(i, it, itrc), itrc
             end do
          end do
       end do
       close(io_syn) 

    end if

    return 
  end subroutine output_results
  
end module mcmc_out
