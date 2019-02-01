!=======================================================================
!   RF_INV: 
!   Trans-dimensional inversion of receiver functions
!   Copyright (C) 2019 Takeshi Akuhara
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!
!   Contact information
!
!   Email  : akuhara @ eri. u-tokyo. ac. jp 
!   Address: Earthquake Research Institute, The Univesity of Tokyo
!           1-1-1, Yayoi, Bunkyo-ku, Tokyo 113-0032, Japan
!
!=======================================================================

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
    integer :: nk_sum(k_max), ierr, ik, itrc, it, i, iv, iz
    integer :: namp_sum(nbin_amp, nsmp, ntrc), nmod_sum
    integer :: nprop_sum(ntype), naccept_sum(ntype)
    integer :: nvsz_sum(nbin_z, nbin_vs), nvpz_sum(nbin_z, nbin_vp)
    integer :: nsig_sum(nbin_sig, ntrc), nz_sum(nbin_z)
    character(clen_max) :: out_file
    
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
    call mpi_reduce(nvpz, nvpz_sum, nbin_z * nbin_vp, &
         & MPI_INTEGER4, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(nvsz, nvsz_sum, nbin_z * nbin_vs, &
         & MPI_INTEGER4, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(nz, nz_sum, nbin_z, &
         & MPI_INTEGER4, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(nsig, nsig_sum, nbin_sig * ntrc, &
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
          write(*,*)"# of perturbing noise sigma proposal:", &
               & naccept_sum(4), "/", nprop_sum(4)
          write(*,*)"# of perturbing Vs proposal:", &
               & naccept_sum(5), "/", nprop_sum(5)
          if (vp_mode == 1) then
             write(*,*)"# of perturbing Vp proposal:", &
                  & naccept_sum(6), "/", nprop_sum(6)
          end if
             write(*,*)
       end if
       
       ! # of layer interfaces
       out_file = trim(out_dir) // "/" // 'num_interface.ppd'
       open(io_nk, file = out_file, status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0,*)"ERROR: cannot create", trim(out_file)
          call mpi_finalize(ierr)
          stop
       end if
       do ik = 1, k_max - 1
          write(io_nk,*) ik, dble(nk_sum(ik)) / dble(nmod_sum)
       end do
       close(io_nk)


       ! Synthetic RFs
       out_file = trim(out_dir) // "/" // 'syn_trace.ppd'
       open(io_syn, file = out_file, status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0,*)"ERROR: cannot create ", trim(out_file)
          call mpi_finalize(ierr)
          stop
       end if
       do itrc = 1, ntrc
          do it = 1, nsmp
             do i = 1, nbin_amp
                write(io_syn,'(3F10.5,I6)') (it-1) * delta + t_start, &
                     & amp_min + (i - 0.5d0) * dbin_amp, &
                     & dble(namp_sum(i, it, itrc)) / dble(nmod_sum), &
                     & itrc
             end do
          end do
       end do
       close(io_syn) 

       ! Interface depth
       out_file = trim(out_dir) // "/" // 'interface_depth.ppd'
       open(io_z, file = out_file, status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0,*)"ERROR: cannot create ", trim(out_file)
          call mpi_finalize(ierr)
          stop
       end if
       do i = 1, nbin_z
          write(io_z,*) &
               & (i - 0.5d0) * dbin_z + z_min, dble(nz_sum(i)) / dble(nmod_sum)
       end do
       close(io_z)

       ! Noise sigma
       out_file = trim(out_dir) // "/" // 'sigma.ppd'
       open(io_sig, file = out_file, status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0,*)"ERROR: cannot create ", trim(out_file)
          call mpi_finalize(ierr)
          stop
       end if
       do itrc = 1, ntrc
          do i = 1, nbin_sig
             write(io_sig,*) (i - 0.5d0) * dbin_sig, &
                  & dble(nsig_sum(i, itrc)) / dble(nmod_sum), &
                  & itrc
          end do
       end do
       close(io_sig)

       ! Vs profile
       out_file = trim(out_dir) // "/" // 'vs_z.ppd'
       open(io_vsz, file = out_file, status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0,*)"ERROR: cannot create ", trim(out_file)
          call mpi_finalize(ierr)
          stop
       end if
       do iv = 1, nbin_vs
          do iz = 1, nbin_z
             write(io_vsz,'(3F10.5)') &
                  & (iv - 0.5d0) * dbin_vs + vs_min, &
                  & (iz - 0.5d0) * dbin_z, &
                  & dble(nvsz_sum(iz, iv)) / dble(nmod_sum)
          end do
       end do
       close(io_vsz)

       ! Vp profile
       out_file = trim(out_dir) // "/" // 'vp_z.ppd'
       open(io_vpz, file = out_file, status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0,*)"ERROR: cannot create ", trim(out_file)
          call mpi_finalize(ierr)
          stop
       end if
       do iv = 1, nbin_vp
          do iz = 1, nbin_z
             write(io_vpz,'(3F10.5)') &
                  & (iv - 0.5d0) * dbin_vp + vp_min, &
                  & (iz - 0.5d0) * dbin_z, &
                  & dble(nvpz_sum(iz, iv)) / dble(nmod_sum)
          end do
       end do
       close(io_vpz)
    end if

    return 
  end subroutine output_results
  
end module mcmc_out
