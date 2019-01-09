module pt_mcmc
  use params
  use mt19937
  implicit none 
  real(8), allocatable :: log_lklh(:)

contains
  
  !---------------------------------------------------------------------
  
  subroutine init_pt_mcmc(verb)
    implicit none 
    logical, intent(in) :: verb
    integer :: ichain
    
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
    
    allocate(log_lklh(nchains))
    
    return 
  end subroutine init_pt_mcmc
  
  !---------------------------------------------------------------------

  subroutine pt_control(verb)
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
    n_all = (nproc - 1) * nchains
    n_tot_iter = nburn + niter
    
    do it = 1, n_tot_iter
       ! Within-chain step for all chains
       if (rank > 0) then
          do ichain = 1, nchains
             temp = temps(ichain)
             call mcmc(ichain, temp, e)
             log_lklh(ichain) = e
          end do
          
          ! Be informed which chains are swapped
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
             e1    = log_lklh(ichain1)
             e2    = log_lklh(ichain2)
             call judge_pt(temp1, temp2, e1, e2, yn)
             if (yn) then
                temps(ichain2) = temp1
                temps(ichain1) = temp2
             end if
          else if (rank1 == rank) then
             ! Receive the other chain's status and judge
             call mpi_recv(rpack, 2, MPI_REAL8, rank2, 2018, &
                  & MPI_COMM_WORLD, status, ierr)
             temp1 = temps(ichain1)
             temp2 = rpack(1)
             e1 = log_lklh(ichain1)
             e2 = rpack(2)
             call judge_pt(temp1, temp2, e1, e2, yn)
             if (yn) then
                temps(ichain1) = temp2
                rpack(1) = temp1
             end if
             call mpi_send(rpack, 1, MPI_REAL8, rank2, 1988, &
                  & MPI_COMM_WORLD, status, ierr)
          else
             ! Send status to the other chain and receive result
             rpack(1) = temps(ichain2)
             rpack(2) = log_lklh(ichain2)
             call mpi_send(rpack, 2, MPI_REAL8, rank1, 2018, &
                  & MPI_COMM_WORLD, status, ierr)
             call mpi_recv(rpack, 1, MPI_REAL8, rank1, 1988, &
                  & MPI_COMM_WORLD, status, ierr)
             temps(ichain2) = rpack(1)
          end if
          
       else
          ! Job of zero-rank processor
          ! - determine chain pair and inform it to all processors
          itarget1 = int(grnd() * n_all)
          do 
             itarget2 = int(grnd() * n_all)
             if (itarget2 /= itarget1) exit 
          end do
          rank1 = int(itarget1 / nchains) + 1
          rank2 = int(itarget2 / nchains) + 1
          ichain1 = mod(itarget1, nchains) + 1
          ichain2 = mod(itarget2, nchains) + 1
          
          ipack(1) = rank1
          ipack(2) = rank2
          ipack(3) = ichain1
          ipack(4) = ichain2
          call mpi_bcast(ipack, 4, MPI_INTEGER4, 0, &
               & MPI_COMM_WORLD, ierr)
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

  subroutine mcmc(ichain, temp, out_log_lklh)
    use model
    use likelihood
    implicit none 
    integer, intent(in) :: ichain
    real(8), intent(in) :: temp
    real(8), intent(out) :: out_log_lklh
    real(8) :: prop_dvp(k_max), prop_dvs(k_max), prop_z(k_max)
    real(8) :: prop_sig(ntrc), prop_log_lklh
    integer :: prop_k
    integer :: itype, itarget, ilay
    logical :: null_flag, yn
    
    prop_k = k(ichain)
    prop_dvp(1:k_max) = dvp(1:k_max, ichain)
    prop_dvs(1:k_max) = dvs(1:k_max, ichain)
    prop_z(1:k_max)   = z(1:k_max, ichain)
    prop_sig(1:ntrc)  = sig(1:ntrc, ichain)

    null_flag = .false.

    ! Select proposal type
    itype = int(grnd() * ntype)

    if (itype == 0) then
       ! Birth proposal
       prop_k = prop_k + 1
       if (prop_k < k_max) then
          prop_dvp(prop_k) = dvp_min + grnd() * (dvp_max - dvp_min)
          prop_dvs(prop_k) = dvs_min + grnd() * (dvs_max - dvs_min)
          prop_z(prop_k)   = z_min   + grnd() * (z_max   - z_min)
       else 
          null_flag = .true.
       end if
    else if (itype == 1) then
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
    else if (itype == 2) then
       ! Move interface
       itarget = int(grnd() * prop_k) + 1
       prop_z(itarget) = prop_z(itarget) + gauss() * dev_z
       if (prop_z(itarget) < z_min .or. &
            & prop_z(itarget) > z_max) then
          null_flag = .true.
       end if
    else if (itype == 3) then
       ! Perturb dVp
       itarget = int(grnd() * (prop_k + 1)) + 1
       prop_dvp(itarget) = prop_dvp(itarget) + gauss() * dev_dvp
       if (prop_dvp(itarget) < dvp_min .or. &
            & prop_dvp(itarget) > dvp_max) then
          null_flag = .true.
       end if
    else if (itype == 4) then
       ! Perturb dVs
       itarget = int(grnd() * (prop_k + 1)) + 1
       prop_dvs(itarget) = prop_dvs(itarget) + gauss() * dev_dvs
       if (prop_dvs(itarget) < dvs_min .or. &
            & prop_dvs(itarget) > dvs_max) then
          null_flag = .true.
       end if
    else 
       ! Perturb noise sigma
       itarget = int(grnd() * ntrc) + 1
       prop_sig(itarget) = prop_sig(itarget) + gauss() * dev_sig
       if (prop_sig(itarget) < sig_min .or. &
            & prop_sig(itarget) > sig_max) then
          null_flag = .true.
       end if
    end if
    
    
    ! evaluate proposed model
    if (.not. null_flag) then
       call calc_log_lklh(ichain, out_log_lklh)
       call judge_mcmc(temp, log_lklh(ichain), out_log_lklh, yn)
       if (yn) then
          log_lklh(ichain) = out_log_lklh
          k(ichain) = prop_k
          dvp(1:k_max, ichain) = prop_dvp(1:k_max)
          dvs(1:k_max, ichain) = prop_dvs(1:k_max)
          z(1:k_max, ichain)   = prop_z(1:k_max)
          sig(1:ntrc, ichain)  = prop_sig(1:ntrc)
       else
          out_log_lklh = log_lklh(ichain)
       end if
    else
       out_log_lklh = log_lklh(ichain)
    end if

    
    return 
  end subroutine mcmc

  !---------------------------------------------------------------------
  subroutine judge_mcmc(temp, log_lklh1, log_lklh2, yn)
    use mt19937
    implicit none
    real(8), intent(in) :: temp, log_lklh1, log_lklh2
    logical, intent(out) :: yn
    real(8) :: del_s
    
    yn = .false.
    del_s = (log_lklh2 - log_lklh1) / temp
    
    if (log(grnd()) <= del_s) then
       yn = .true.
    end if
    
    return 
  end subroutine judge_mcmc
  
  !---------------------------------------------------------------------

  real(8) function gauss()
    use mt19937
    implicit none 
    real(8), parameter :: pi2 = 2.d0 * 3.1415926535897931d0
    real(8) :: v1, v2
    
    v1 = grnd() 
    v2 = grnd()
    gauss = sqrt(-2.d0 * log(v1)) * cos(pi2 * v2)
    
    return 
  end function gauss
  
  !---------------------------------------------------------------------


end module pt_mcmc
