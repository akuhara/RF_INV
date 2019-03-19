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

module pt_mcmc
  use params
  implicit none 
  integer :: nmod 
  integer, allocatable, public :: nk(:), nsig(:,:), namp(:,:,:), nz(:)
  integer, allocatable, public :: nprop(:), naccept(:)
  integer, allocatable, public :: nvpz(:,:), nvsz(:,:)
  real(8), allocatable, public :: likelihood_hist(:)
  real(8), public :: dbin_vp, dbin_vs, dbin_z, dbin_amp, dbin_sig

  integer, public :: ntype, itype_birth, itype_death
  integer, public :: itype_dvs, itype_dvp, itype_sig, itype_z
  character(clen_max), allocatable, public :: prop_label(:)

  public init_pt_mcmc
  private judge_mcmc, judge_pt, mcmc
contains
  
  !---------------------------------------------------------------------
  subroutine mcmc(iter, ichain, temp)
    use mt19937
    use model
    use likelihood
    use forward
    use prior
    use math
    implicit none 
    integer, intent(in) :: iter, ichain
    real(8), intent(in) :: temp
    real(8) :: prop_log_likelihood, prop_rft(nfft, ntrc)
    real(8) :: log_prior12
    real(8) :: prop_dvp(k_max), prop_dvs(k_max), prop_z(k_max)
    real(8) :: prop_sig(ntrc), tmpz
    integer :: prop_k
    integer :: itype, itarget, ilay, ibin, iz1, iz2, itrc, iz
    integer :: ivp, ivs, it
    logical :: null_flag, yn, fwd_flag, is_valid
    integer :: nlay
    real(8) :: alpha(nlay_max), beta(nlay_max)
    real(8) :: h(nlay_max), rho(nlay_max)

    
    prop_k = k(ichain)
    prop_dvp(1:k_max)   = dvp(1:k_max, ichain)
    prop_dvs(1:k_max)   = dvs(1:k_max, ichain)
    prop_z(1:k_max - 1) = z(1:k_max - 1, ichain)
    prop_sig(1:ntrc)    = sig(1:ntrc, ichain)

    null_flag = .false.

    ! Select proposal type
    itype = int(grnd() * ntype) + 1

    ! Make Proposal

    if (itype == itype_birth) then
       ! Birth proposal
       prop_k = prop_k + 1
       if (prop_k < k_max) then
          prop_dvp(prop_k) = gauss() * dvp_prior
          prop_dvs(prop_k) = gauss() * dvs_prior
          prop_z(prop_k)   = z_min   + grnd() * (z_max   - z_min)
       else 
          null_flag = .true.
       end if
       
    else if (itype == itype_death) then
       ! Death proposal
       prop_k = prop_k - 1
       if (prop_k >= k_min) then
          itarget = int(grnd() * (prop_k + 1)) + 1
          do ilay = itarget, prop_k
             prop_dvp(ilay) = dvp(ilay + 1, ichain)
             prop_dvs(ilay) = dvs(ilay + 1, ichain)
             prop_z(ilay)   = z(ilay + 1, ichain)
          end do
          prop_dvp(prop_k + 1) = 0.d0
          prop_dvs(prop_k + 1) = 0.d0
          prop_z(prop_k + 1)   = 0.d0
       else
          null_flag = .true.
       end if
       
    else if (itype == itype_z) then
       ! Move interface
       itarget = int(grnd() * prop_k) + 1
       prop_z(itarget) = prop_z(itarget) + gauss() * dev_z
       if (prop_z(itarget) < z_min .or. &
            & prop_z(itarget) > z_max) then
          null_flag = .true.
       end if
       
    else if (itype == itype_dvs) then
       ! Perturb dVs
       itarget = int(grnd() * (prop_k + 1)) + 1
       if (itarget == prop_k + 1) itarget = k_max
       prop_dvs(itarget) = prop_dvs(itarget) + gauss() * dev_dvs
       log_prior12 = &
            & log_prior_ratio(prop_dvs(itarget), &
            & dvs(itarget, ichain), dvs_prior)
       
    else if (itype == itype_dvp) then
       ! Perturb dVp
       itarget = int(grnd() * (prop_k + 1)) + 1
       if (itarget == prop_k + 1) itarget = k_max
       prop_dvp(itarget) = prop_dvp(itarget) + gauss() * dev_dvp
       log_prior12 = &
            & log_prior_ratio(prop_dvp(itarget), &
            & dvp(itarget, ichain), dvp_prior)

    else if (itype == itype_sig) then
       ! Perturb noise sigma
       itarget = int(grnd() * ntrc) + 1
       prop_sig(itarget) = prop_sig(itarget) + gauss() * dev_sig
       if (prop_sig(itarget) < sig_min .or. &
            & prop_sig(itarget) > sig_max) then
          null_flag = .true.
       end if
    end if
    
    ! Make velocity model
    if (.not. null_flag) then
       call format_model(prop_k, prop_z, prop_dvp, prop_dvs, &
            & nlay, alpha, beta, rho, h, is_valid)
    end if
    if (.not. is_valid) then ! check velocity range
       null_flag = .true.
    end if
    
    ! evaluate proposed model
    if (.not. null_flag) then
       if (itype /= 4) then
          fwd_flag = .true.
       else
          fwd_flag = .false.
       end if
       call calc_likelihood(ichain, fwd_flag, &
            & prop_k, prop_z, prop_dvp, prop_dvs, &
            & prop_sig, prop_log_likelihood, prop_rft)
       call judge_mcmc(temp, log_likelihood(ichain), &
            & prop_log_likelihood, log_prior12, yn)
       if (yn) then
          log_likelihood(ichain) = prop_log_likelihood
          k(ichain)              = prop_k
          dvp(1:k_max, ichain)   = prop_dvp(1:k_max)
          dvs(1:k_max, ichain)   = prop_dvs(1:k_max)
          z(1:k_max - 1, ichain) = prop_z(1:k_max - 1)
          sig(1:ntrc, ichain)    = prop_sig(1:ntrc)
          rft(1:nfft, 1:ntrc, ichain) = prop_rft(:,:)
       end if
    end if

    
    ! count proposal
    if (temp <= 1.d0 + 1.0e-6) then
       nprop(itype) = nprop(itype) + 1
       if (yn) naccept(itype) = naccept(itype) + 1

       likelihood_hist(iter) = &
            likelihood_hist(iter) + log_likelihood(ichain)
       
    end if

    ! record sampled model
    if (temp <= 1.d0 + 1.0e-6 .and. iter > nburn .and. &
         & mod(iter, ncorr) == 0) then
       
       nmod = nmod + 1

       
       ! # of layer interfaces
       nk(k(ichain)) = nk(k(ichain)) + 1
       
       ! Noise sigma
       do itrc = 1, ntrc
          ibin = int((sig(itrc, ichain) - sig_min) / dbin_sig) + 1
          nsig(ibin, itrc) = nsig(ibin, itrc) + 1
       end do

       ! Interface depth
       do ilay = 1, k(ichain) - 1
          ibin = int((z(ilay, ichain) - z_min) / dbin_z) + 1
          nz(ibin) = nz(ibin) + 1
       end do


       ! V-z profile
       call format_model(k(ichain), z(1:k_max-1, ichain), &
            & dvp(1:k_max, ichain), dvs(1:k_max, ichain), &
            & nlay, alpha, beta, rho, h, is_valid)
       tmpz = 0.d0
       do ilay = 1, nlay
          iz1 = int(tmpz / dbin_z) + 1
          if (ilay < nlay) then
             iz2 = int((tmpz + h(ilay)) / dbin_z) + 1
          else
             iz2 = nbin_z + 1
          end if
          ivp = int((alpha(ilay) - vp_min) / dbin_vp) + 1
          ivs = int((beta(ilay) - vs_min) / dbin_vs) + 1
          ivs = max(1, ivs) ! for ocean layer where beta(1) < vs_min
          do iz = iz1, iz2 - 1
             nvpz(iz, ivp) = nvpz(iz, ivp) + 1
             nvsz(iz, ivs) = nvsz(iz, ivs) + 1
          end do
          tmpz = tmpz + h(ilay)
       end do

       ! RF trace
       do itrc = 1, ntrc
          do it = 1, nsmp
             ibin = int((rft(it, itrc, ichain) - amp_min) &
                  & / dbin_amp) + 1
             if (ibin < 1) then
                write(0,*)"Warning: RF amp. out of range"
                ibin = 1
             else if (ibin > nbin_amp) then
                write(0,*)"Warning: RF amp. out of range"
                ibin = nbin_amp
             end if
             namp(ibin, it, itrc) = namp(ibin, it, itrc) + 1
          end do
       end do
    end if

    
    return 
  end subroutine mcmc



  !---------------------------------------------------------------------
  
  subroutine init_pt_mcmc(verb)
    use model
    use likelihood
    implicit none 
    logical, intent(in) :: verb
    real(8) :: prop_rft(nfft, ntrc)
    integer :: ichain, ierr, itype

    !allocate
    if (verb) then
       write(*,*)
       write(*,*) "--- Allocate memory for recording models ---"
    end if
    
    allocate(nk(k_max), nz(nbin_z), nsig(nbin_sig, ntrc))
    allocate(namp(nbin_amp, nsmp, ntrc))
    allocate(nvpz(nbin_z, nbin_vp), nvsz(nbin_z, nbin_vs))
    allocate(nprop(ntype), naccept(ntype))
    allocate(likelihood_hist(niter + nburn))
    
    nk = 0
    nz = 0
    nsig = 0
    namp = 0
    nvpz = 0
    nvsz = 0
    nmod = 0
    naccept = 0
    nprop = 0
    likelihood_hist = 0.d0

    dbin_sig = (sig_max - sig_min) / nbin_sig
    dbin_amp = (amp_max - amp_min) / nbin_amp
    dbin_vp = (vp_max - vp_min) / nbin_vp
    dbin_vs = (vs_max - vs_min) / nbin_vs
    dbin_z = (z_max - 0.d0) / nbin_z
    


    if (verb) then
       write(*,*) "  Done"
    end if
    
    
    if (verb) then
       write(*,*)
       write(*,*) "--- Initialize temperatures ---"
    end if

    allocate(temps(nchains))

    ! non-tempered chain
    temps(1:ncool) = 1.d0

    ! tempered chain
    do ichain = ncool+1, nchains
       temps(ichain) = exp(grnd() * log(t_high))
    end do
  
    if (verb) then
       do ichain = 1, nchains
          write(*,*) temps(ichain), &
               & "Note: This is a sample value - not used actually."
       end do
    end if
    
    ! Proposal type
    if (verb) then
       write(*,*)
       write(*,*) "Now checking proposal type"
    end if
    ntype = 4
    itype_birth = 1
    itype_death = 2
    itype_z     = 3
    itype_dvs   = 4
    if (vp_mode == 1) then
       ntype = ntype + 1
       itype_dvp = ntype
       if (verb) then
          write(*,*)"dVp is solved"
       end if
    else if (vp_mode == 0) then
       itype_dvp = -1
       if (verb) then 
          write(*,*)"dVp is fixed at 0"
       end if
    else
       if (verb) then
          write(*,*)"ERRPR: vp_mode sholud be 0 or 1"
          write(*,*)  " 0: fixed, 1: solved"
          call mpi_finalize(ierr)
          stop
       end if
    end if
    if (sig_mode == 1) then
       ntype = ntype + 1
       itype_sig = ntype
    else
       itype_sig = -1
    end if

    ! Make proposal label for output summary 
    allocate(prop_label(ntype))
    do itype = 1, ntype
       if (itype == itype_birth) then
          prop_label(itype) = "Birth proposal"
       else if (itype == itype_death) then
          prop_label(itype) = "Death proposal"
       else if (itype == itype_z) then
          prop_label(itype) = "Moving interface depth proposal"
       else if (itype == itype_dvs) then
          prop_label(itype) = "Perturbing dVs proposal"
       else if (itype == itype_dvp) then
          prop_label(itype) = "Perturbing dVs proposal"
       else if (itype == itype_sig) then
          prop_label(itype) = "Perturbing sigma proposal"
       end if
    end do
    

    return 
  end subroutine init_pt_mcmc
  
  !---------------------------------------------------------------------

  subroutine pt_control(verb)
    use likelihood
    use mt19937
    implicit none 
    include "mpif.h"
    logical, intent(in) :: verb
    integer :: n_all, n_tot_iter, nproc, ierr, rank
    integer :: it, ichain, ipack(4), rank1, rank2
    integer :: ichain1, ichain2, status(MPI_STATUS_SIZE)
    integer :: itarget1, itarget2
    real(8) :: temp, e, temp1, temp2, e1, e2, rpack(2)
    logical :: yn

    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank,  ierr)
    

    n_all = nproc * nchains
    n_tot_iter = nburn + niter
    
    do it = 1, n_tot_iter
       if (rank == 1 .and. mod(it, ncorr) == 0) then
          write(*,*)"Iteration #:", it, "/", n_tot_iter
       end if
       ! Within-chain step for all chains
       do ichain = 1, nchains
          temp = temps(ichain)
          call mcmc(it, ichain, temp)
       end do
       
       if (nchains < 2) cycle ! single-chain MCMC

       ! determine chain pair (extra job for rank=0 chain)
       if (rank == 0) then
          itarget1 = int(grnd() * n_all)
          do 
             itarget2 = int(grnd() * n_all)
             if (itarget2 /= itarget1) exit 
          end do

          rank1 = int(itarget1 / nchains)
          rank2 = int(itarget2 / nchains)
          ichain1 = mod(itarget1, nchains) + 1
          ichain2 = mod(itarget2, nchains) + 1
          
          ipack(1) = rank1
          ipack(2) = rank2
          ipack(3) = ichain1
          ipack(4) = ichain2  
       end if
       call mpi_bcast(ipack, 4, MPI_INTEGER4, 0, &
            & MPI_COMM_WORLD, ierr)
       rank1   = ipack(1)
       rank2   = ipack(2)
       ichain1 = ipack(3)
       ichain2 = ipack(4)
       
       if (rank1 == rank .and. rank2 == rank) then
          ! Swap within the single processor
          temp1 = temps(ichain1)
          temp2 = temps(ichain2)
          e1    = log_likelihood(ichain1)
          e2    = log_likelihood(ichain2)
          call judge_pt(temp1, temp2, e1, e2, yn)
          if (yn) then
             temps(ichain2) = temp1
             temps(ichain1) = temp2
          end if
          !if (yn .and. temp1 < 1.d0 + 1.0e-6) then
          !   write(*,*)"Successful jump of non-tempered chain to ", temp2
          !else if (yn .and. temp2 < 1.d0 + 1.0e-6) then
          !   write(*,*)"Successful jump of non-tempered chain to ", temp1
          !end if
          
       else if (rank1 == rank) then
          ! Receive the other chain's status and judge
          call mpi_recv(rpack, 2, MPI_REAL8, rank2, 2018, &
               & MPI_COMM_WORLD, status, ierr)
          temp1 = temps(ichain1)
          temp2 = rpack(1)
          e1 = log_likelihood(ichain1)
          e2 = rpack(2)
          call judge_pt(temp1, temp2, e1, e2, yn)
          if (yn) then
             temps(ichain1) = temp2
             rpack(1) = temp1
          end if
          call mpi_send(rpack, 1, MPI_REAL8, rank2, 1988, &
               & MPI_COMM_WORLD, status, ierr)
          !if (yn .and. temp1 < 1.d0 + 1.0e-6) then
          !   write(*,*)"Successful jump of non-tempered chain to ", temp2
          !else if (yn .and. temp2 < 1.d0 + 1.0e-6) then
          !   write(*,*)"Successful jump of non-tempered chain to ", temp1
          !end if
       else if (rank2 == rank) then
          ! Send status to the other chain and receive result
          rpack(1) = temps(ichain2)
          rpack(2) = log_likelihood(ichain2)
          call mpi_send(rpack, 2, MPI_REAL8, rank1, 2018, &
               & MPI_COMM_WORLD, status, ierr)
          call mpi_recv(rpack, 1, MPI_REAL8, rank1, 1988, &
               & MPI_COMM_WORLD, status, ierr)
          temps(ichain2) = rpack(1)
       end if
    end do
    
       

  end subroutine pt_control
  
  !---------------------------------------------------------------------

  subroutine judge_pt(temp1, temp2, l1, l2, yn)
    use mt19937
    implicit none 
    real(kind(0d0)), intent(in) :: temp1, temp2, l1, l2
    logical, intent(out) :: yn
    real(kind(0d0)) :: del_s
    
    del_s = (l2 - l1) * (1.d0 / temp1 - 1.d0 / temp2)
    yn = .false.
    
    if(log(grnd()) <= del_s) then
       yn = .true.
    end if
    
    return 
  end subroutine judge_pt

  !---------------------------------------------------------------------

 
  subroutine judge_mcmc(temp, log_lklh1, log_lklh2, log_prior12, yn)
    use mt19937
    implicit none
    real(8), intent(in) :: temp, log_lklh1, log_lklh2, log_prior12
    logical, intent(out) :: yn
    real(8) :: del_s
    
    yn = .false.
    del_s = (log_lklh2 - log_lklh1) / temp + log_prior12
    
    if (log(grnd()) <= del_s) then
       yn = .true.
    end if
    
    return 
  end subroutine judge_mcmc
  
  !---------------------------------------------------------------------

  


end module pt_mcmc
