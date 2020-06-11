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
  
  subroutine output_results(nproc, rank, verb)
    implicit none 
    include "mpif.h"
    integer, intent(in) :: nproc, rank
    logical, intent(in) :: verb
    integer :: nk_sum(k_max), ierr, ik, itrc, it, i, iv, iz, imod
    integer :: namp_sum(nbin_amp, nsmp, ntrc), nmod_sum
    integer :: nprop_sum(ntype), naccept_sum(ntype)
    integer :: nvsz_sum(nbin_z, nbin_vs), nvpz_sum(nbin_z, nbin_vp)
    integer :: nvpvsz_sum(nbin_z, nbin_vpvs)
    integer :: nsig_sum(nbin_sig, ntrc), nz_sum(nbin_z)
    real(8) :: likelihood_hist_av(nburn + niter)
    real(8) :: vp_mean_sum(nbin_z), vs_mean_sum(nbin_z), vpvs_mean_sum(nbin_z)
    real(8), allocatable :: vp_model_sum(:,:), vs_model_sum(:,:)
    real(8), allocatable :: all_likelihood_sum(:)
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
    call mpi_reduce(nvpvsz, nvpvsz_sum, nbin_z * nbin_vpvs, &
         & MPI_INTEGER4, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(nz, nz_sum, nbin_z, &
         & MPI_INTEGER4, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(nsig, nsig_sum, nbin_sig * ntrc, &
         & MPI_INTEGER4, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(likelihood_hist, likelihood_hist_av, &
         & nburn + niter, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(vp_mean, vp_mean_sum, nbin_z, MPI_REAL8, &
         & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(vs_mean, vs_mean_sum, nbin_z, MPI_REAL8, &
         & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(vpvs_mean, vpvs_mean_sum, nbin_z, MPI_REAL8, &
         & MPI_SUM, 0, MPI_COMM_WORLD, ierr)


    
    allocate(all_likelihood_sum(int(niter*nchains*nproc/ncorr)))
    allocate(vp_model_sum(nbin_z, int(niter*nchains*nproc/ncorr)))
    allocate(vs_model_sum(nbin_z, int(niter*nchains*nproc/ncorr)))
    
    
    call mpi_gather(vs_model, nbin_z*int(niter*nchains/ncorr), &
         & MPI_REAL8, vs_model_sum, nbin_z*int(niter*nchains/ncorr), &
         & MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call mpi_gather(vp_model, nbin_z*int(niter*nchains/ncorr), &
         & MPI_REAL8, vp_model_sum, nbin_z*int(niter*nchains/ncorr), &
         & MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    

    if (rank == 0) then
       
       if (verb) then
          write(*,*)"--- Summary ---"
          write(*,*)"# of sampled models:", nmod_sum
          do i = 1, ntype
             write(*,*)"# of ", trim(prop_label(i)), ": ", &
             & naccept_sum(i), "/", nprop_sum(i)
          end do
          write(*,*)
       end if
       
       ! All models
       out_file = trim(out_dir) // "/" // 'all_models'
       open(141, file = out_file, status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0,*)"ERROR: cannot create", trim(out_file)
          call mpi_finalize(ierr)
          stop
       end if
       do imod = 1, niter*nchains*nproc/ncorr
          if (vs_model_sum(1, imod) < -900.d0) then
             cycle
          end if
          write(141,*)""
          do iz = 1, nbin_z
             write(141,*)(iz - 0.5d0) * dbin_z, vp_model_sum(iz, imod), &
                  & vs_model_sum(iz, imod)
          end do
          write(141,*)""
       end do
       
       close(141)

       ! Transition of Likelihood 
       out_file = trim(out_dir) // "/" // 'likelihood'
       open(io_lkhd, file = out_file, status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0,*)"ERROR: cannot create", trim(out_file)
          call mpi_finalize(ierr)
          stop
       end if
       do it = 1, nburn + niter
          write(io_lkhd,*)it, likelihood_hist_av(it) / dble(ncool * nproc)
       end do
       close(io_lkhd)

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
               & (i - 0.5d0) * dbin_z, dble(nz_sum(i)) / dble(nmod_sum)
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
          if (sig_mode(itrc) == 1) then
             do i = 1, nbin_sig
                write(io_sig,*) (i - 0.5d0) * dbin_sig(itrc) &
                     & + sig_min(itrc), &
                     & dble(nsig_sum(i, itrc)) / dble(nmod_sum), &
                     & itrc
             end do
          end if
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
       
       ! Vp/Vs profile
       out_file = trim(out_dir) // "/" // 'vpvs_z.ppd'
       open(io_vpvsz, file = out_file, status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0,*)"ERROR: cannot create ", trim(out_file)
          call mpi_finalize(ierr)
          stop
       end if
       do iv = 1, nbin_vpvs
          do iz = 1, nbin_z
             write(io_vpvsz,'(3F10.5)') &
                  & (iv - 0.5d0) * dbin_vpvs + vpvs_min, &
                  & (iz - 0.5d0) * dbin_z, &
                  & dble(nvpvsz_sum(iz, iv)) / dble(nmod_sum)
          end do
       end do
       close(io_vpvsz)

       ! Vs mean
       out_file = trim(out_dir) // "/" // "vs_z.mean"
       open(io_vs_mean, file = out_file, status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0,*)"ERROR: cannot create ", trim(out_file)
          call mpi_finalize(ierr)
          stop
       end if
       do iz = 1, nbin_z
          write(io_vs_mean,'(2F10.5)')  dble(vs_mean_sum(iz)) / &
               & dble(nmod_sum), (iz - 0.5d0) * dbin_z
       end do
       close(io_vs_mean)
       
       ! Vp mean
       out_file = trim(out_dir) // "/" // "vp_z.mean"
       open(io_vp_mean, file = out_file, status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0,*)"ERROR: cannot create ", trim(out_file)
          call mpi_finalize(ierr)
          stop
       end if
       do iz = 1, nbin_z
          write(io_vp_mean,'(2F10.5)')  dble(vp_mean_sum(iz)) / &
               & dble(nmod_sum), (iz - 0.5d0) * dbin_z
       end do
       close(io_vp_mean)
       
       ! Vp/Vs mean
       out_file = trim(out_dir) // "/" // "vpvs_z.mean"
       open(io_vpvs_mean, file = out_file, status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0,*)"ERROR: cannot create ", trim(out_file)
          call mpi_finalize(ierr)
          stop
       end if
       do iz = 1, nbin_z
          write(io_vpvs_mean,'(2F10.5)')  dble(vpvs_mean_sum(iz)) / &
               & dble(nmod_sum), (iz - 0.5d0) * dbin_z
       end do
       close(io_vpvs_mean)
    end if

    return 
  end subroutine output_results
  
end module mcmc_out
