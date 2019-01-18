module PTmod
  implicit none 
  !        nproc      Integer   : Number of processes used (MPI mode only)
  !        rank       Integer   : Rank of process (MPI mode only)
  !        iseed      Integer   : Seed for random number generator
  !        dir        Character : Base directory for I/O (Some mpi jobs require full path names)
  !        AllTemps   Double    : Temperature list across all chains on all processes 
  !
  Integer                            :: nproc,rank
  Character(len=80)                     dir
  Double precision, allocatable      :: AllTemps(:)
   
end module PTmod
!----------------------------------------------------------------------
!-------------------------------------------------------------------------
!
!     Global data arrays defining model space - used by driver program 
!
!-------------------------------------------------------------------------
  
module modelspacedata
  Real(8), allocatable :: model(:,:) !state of each chain model space
  Real(8), allocatable :: propmodel(:) !proposed model array 
  Double precision, allocatable :: logPPDstore(:) !current logPPD of each chain 
  real(8), allocatable :: rft_store(:,:,:)
  Integer, allocatable :: n_vz(:,:),  n_vpvsz(:,:)
  Integer, allocatable :: n_int(:)
  Integer, allocatable :: n_amp(:,:,:)
  Integer, allocatable :: n_dim(:), n_sig(:,:)
  Integer, allocatable :: n_propose(:), n_accept(:)
  Integer :: ntransprop,ntrans, nmod
  Integer :: iburnMC,icountMcMC
  integer, allocatable :: step_count(:)
end module modelspacedata

!------------------------------------------------------------


program inv_PT_RF
  use PTmod
  use modelspacedata
  use params
  implicit none 
#if defined MPI
  include "mpif.h"
#endif
  Double precision, allocatable      :: Tbins(:)    
  Integer, allocatable               :: Ttot(:)     
  Double precision, allocatable      :: Chaintemp(:)
  Double precision                   :: Tlow
  Real, allocatable                  :: temp_accept_rates(:,:)
  Real, allocatable                  :: mcmc_accept_rates(:)
  Real, allocatable                  :: sjd(:,:)
  Integer, parameter              :: ivzfile = 10, ivzhotfile= 20, iintfile = 30, iampfile=40, ivpvszfile=70
  Integer, parameter              :: idimfile=50, isigfile=60, imodfile=100
  Integer, allocatable :: n_vz_sum(:,:), n_int_sum(:), n_amp_sum(:,:,:)
  Integer, allocatable :: n_dim_sum(:), n_vpvsz_sum(:,:)
  Integer, allocatable :: n_sig_sum(:,:), nmods(:)
  Integer, allocatable :: n_propose_sum(:), n_accept_sum(:)
  integer :: ialg, nbins
  real :: swaprate
  real(8) :: z, v, t, a, h
  integer :: i, j, nerr, nmod_sum, k
  character(100) :: filename, str

  

  
  ialg     = 2            ! Algorithm type: 
  ! ialg=0, Parallel Tempering with T swap between all levels 
  ! ialg=2, Parallel Tempering with T swap only allowed between 
  ! neighbouring temperature levels 
  
  Tlow  = 1.d0            
  nbins = 10   !10       
  dir = './'             
  
  
  call pt(0,0,0,0,0,0.D0,0.D0,0.D0,0,0.0,0,dir,nproc,rank)
  allocate(nmods(nproc))


  call Initialize_pspace()  
  if (ncool < nchains) then
     swaprate = 1.d0     
  else
     swaprate = 0.d0
  end if


  allocate(Chaintemp(nchains)) 
  allocate(Tbins(nbins))           
  allocate(Ttot(nbins))  
  allocate(temp_accept_rates(nbins,nbins))
  allocate(mcmc_accept_rates(nbins))      
  allocate(sjd(0:nbins,nbins))            
  

  call Setuptempladder(nchains,ncool,Tlow,Thigh,Chaintemp) 
  call Setuptempbins(nbins,Tlow,Thigh,nchains,Tbins,Ttot)  
  call PT_diagnostics(0,nbins,Tbins,temp_accept_rates,mcmc_accept_rates,sjd)   
  call PT &
       (1,ialg,nchains,nsteps,iburn,&
       &Chaintemp,Thigh,Tlow,nbins,swaprate,&
       &iseed0,dir,nproc,rank)
  call PT_diagnostics(1,nbins,Tbins,temp_accept_rates,mcmc_accept_rates,sjd)   

  

! output
#if defined MPI
  
  write(*,*)"REDUCE"
  write(*,*)"nmod", nmod, rank
  allocate(n_vz_sum(n_vbin,n_zbin))
  allocate(n_vpvsz_sum(n_vpvsbin,n_zbin))
  allocate(n_int_sum(n_zbin))
  allocate(n_amp_sum(n_smp,n_abin,ntrc))
  allocate(n_dim_sum(kmax-1))
  allocate(n_sig_sum(n_sigbin,ntrc))
  allocate(n_propose_sum(ntype))
  allocate(n_accept_sum(ntype))
  



  call MPI_Reduce(n_vz(1,1),n_vz_sum(1,1),n_vbin*n_zbin,MPI_INTEGER4, &
       & MPI_SUM,0,MPI_COMM_WORLD,nerr)
  call MPI_Reduce(n_vpvsz(1,1),n_vpvsz_sum(1,1),n_vpvsbin*n_zbin,MPI_INTEGER4, &
       & MPI_SUM,0,MPI_COMM_WORLD,nerr)
  call MPI_Reduce(n_int(1),n_int_sum(1),n_zbin,MPI_INTEGER4, &
       & MPI_SUM,0,MPI_COMM_WORLD,nerr)
  call MPI_Reduce(n_amp(1,1,1),n_amp_sum(1,1,1),n_abin*n_smp*ntrc,&
       & MPI_INTEGER4, MPI_SUM,0,MPI_COMM_WORLD,nerr)
  call MPI_Reduce(n_dim(1),n_dim_sum(1),kmax-1,MPI_INTEGER4, &
      & MPI_SUM,0,MPI_COMM_WORLD,nerr)
  call MPI_Reduce(n_sig(1,1),n_sig_sum(1,1),n_sigbin*ntrc,MPI_INTEGER4, &
       & MPI_SUM,0,MPI_COMM_WORLD,nerr)
  call MPI_Reduce(nmod,nmod_sum,1,MPI_INTEGER4, &
       & MPI_SUM,0,MPI_COMM_WORLD,nerr)
  !call MPI_Gather(nmod,1,MPI_INTEGER4,nmods,1,MPI_INTEGER4,0,&
  !     & MPI_COMM_WORLD, nerr)
  call MPI_Reduce(n_accept(1),n_accept_sum(1),ntype,MPI_INTEGER4,&
       & MPI_SUM,0,MPI_COMM_WORLD,nerr)
  call MPI_Reduce(n_propose(1),n_propose_sum(1),ntype,MPI_INTEGER4,&
       & MPI_SUM,0,MPI_COMM_WORLD,nerr)
  


  if (rank == 0) then
     write(*,*)
     do i = 1, ntype
        if (i == itype_vs) then
           str = "Vs"
        else if (i == itype_z) then
           str = "Depth"
        else if (i == itype_vpvs) then
           str = "Vp/Vs"
        else if (i == itype_sig) then
           str = "sigma"
        else if (i == itype_birth) then
           str = "birth"
        else if (i == itype_death) then
           str = "death"
        else
           str = "ERROR"
        end if
        write(*,*) "proposal type: ", trim(str), " ", &
             & n_accept_sum(i), "/", n_propose_sum(i)
     end do
     write(*,*)

     ! output 
     !** V-Z
     open(ivzfile,file="vs_z.dat",status='replace')
     do i = 1, n_zbin
        z = (i-1+0.5) * del_zbin
        do j = 1, n_vbin
           v = vbinmin + (j-1+0.5) * del_vbin
          write(ivzfile,'(3F13.6)')v, z, &
                & dble(n_vz_sum(j,i)) / dble(nmod_sum)
        end do
     end do
     close(ivzfile) 

     !** vpvs-Z
     open(ivpvszfile,file="vpvs_z.dat",status='replace')
     do i = 1, n_zbin
        z = (i-1+0.5) * del_zbin
        do j = 1, n_vpvsbin
           v = vpvsmin + (j-1+0.5) * del_vpvsbin
          write(ivpvszfile,'(3F13.6)')v, z, &
                & dble(n_vpvsz_sum(j,i)) / dble(nmod_sum)
        end do
     end do
     close(ivpvszfile) 
     

     !** Layer interface depths
     open(iintfile,file="interface_depth.dat",status="replace")
     do i = 1, n_zbin
        z = (i-1+0.5) * del_zbin + sdep
        write(iintfile,'(2F13.6)') z, dble(n_int_sum(i)) / dble(nmod_sum)
     end do
     close(iintfile)
     
     !** Synthetic RFs
     do k = 1, ntrc
        write(filename,'(A6,I2.2,A4)')"syn_rf", k, ".dat"
        open(iampfile+k,file=filename,status="replace")
        do i = 1, n_abin
           a = abinmin + (i-1 + 0.5) * del_abin 
           do j = 1, n_smp
              t = t_win_start + (j-1)*delta
              write(iampfile+k,'(3F13.6)')t, a, dble(n_amp_sum(j,i,k)) / dble(nmod_sum)
           end do
        end do
        close(iampfile+k)
     end do
  
     !** # of layers
     open(idimfile,file="dimension.dat",status="replace")
     do i = 1, kmax-1
        write(idimfile,'(2F13.6)')dble(i), dble(n_dim_sum(i)) / dble(nmod_sum)
     end do
     close(idimfile)
     
     !** sigma
     do k = 1, ntrc
        write(filename,'(A4,I2.2,A4)')"stdv", k, ".dat"
        open(isigfile+k,file=filename,status="replace")
        do i = 1, n_sigbin
           write(isigfile+k,'(2F13.6)')dble(i-1+0.5)*del_sigbin + sigmin, dble(n_sig_sum(i,k)) / dble(nmod_sum)
        end do
        close(isigfile+k)
     end do
  end if
#endif

  call MPI_Finalize(nerr)
  !call pt (99,0,0,0,0,0.D0,0.D0,0.D0,0,0.0,0,dir,nproc,rank)
  write(*,*)"FINISH", rank
  

  stop
end program inv_PT_RF
!----------------------------------------------------------------------

!---------------------------------------------------------------------
subroutine Initialize_pspace()
  use PTmod
  use modelspacedata
  use params
  use mt19937
  implicit none 
  integer :: i, j, k, itrc
  real(8) :: rand
  

  write(*,*)"Initialize"

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Read parameter file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call read_param("params.in")
  call calc_param()
  
  allocate(model(npara,nchains))         ! Allocate model space arrays
  allocate(propmodel(npara))
  allocate(logPPDstore(nchains))
  allocate(rft_store(nfft,ntrc,nchains))
  allocate(n_vz(n_vbin,n_zbin))
  allocate(n_vpvsz(n_vpvsbin,n_zbin))
  allocate(n_int(n_zbin))
  allocate(n_amp(n_smp,n_abin,ntrc))
  allocate(n_dim(kmax-1))
  allocate(n_sig(n_sigbin,ntrc))
  allocate(n_propose(ntype), n_accept(ntype))
  allocate(step_count(nchains))
  write(*,*)"Finish allocating memory"

  ntransprop = 0                       ! Initialize variables for acceptance probability of jumps
  ntrans = 0                           ! Initialize variables for acceptance probability of jumps
  icountMcMC = 0                       ! Initialize some variables for routine McMC
  iburnMC = iburn*nchains 
  model = 0.d0
  propmodel = 0.d0                     ! initialize model perturbations
  n_int = 0
  n_vz = 0
  n_vpvsz= 0
  n_amp    = 0
  n_dim = 0
  n_sig = 0
  nmod = 0
  n_accept = 0
  n_propose = 0
  step_count = 0
  iseed = -(iseed + 1000 * rank * rank)
  call sgrnd(iseed)
  

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Read reference velocity file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call read_vel_file()
  write(*,*)"Finish reading reference velocity model"
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Read data file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call read_data_file()
  write(*,*)"Finish reading observed RF data"
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Calculate R inverse
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call calc_rinv2()

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! generate initial model in random manner 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do i = 1, nchains
     ! set dimension
     k = kmin + int(grnd() * del_k)
     model(1,i) = dble(k)

     ! set interface-depth 
     do j = 1, k
        model(j+1,i) = zmin + grnd() * del_z
     end do

     ! set velocity
     do j = 1, k 
        model(j+kmax+1,i) = vmin + grnd()*del_v
     end do
     model(2*kmax+1,i) = vmin + grnd()*del_v ! bottom half space

     ! set Vp/Vs ratio
     do j = 1, k
        model(j+2*kmax+1,i) = vpvsmin + grnd()*del_vpvs
     end do
     model(3*kmax+1,i) = vpvsmin + grnd()*del_vpvs ! bottom half space
     
     ! set noise sigma
     do itrc = 1, ntrc
        model(3*kmax+1+itrc,i) = sigmin + grnd()*del_sig
     end do
     
     ! sediment
     if (sed_flag == 1) then
        model(3*kmax+ntrc+2,i) = sed_vsmin + grnd()*del_sedvs
        model(3*kmax+ntrc+3,i) = sed_vpmin + grnd()*del_sedvp
        model(3*kmax+ntrc+4,i) = sed_hmin + grnd()*del_sedh
     end if
     
  end do
  write(*,*)"Finish generating initial model"




  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! calculate log PPD
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  do i = 1, nchains
     k = int(model(1,i))
     call CalclogPPD(k, model(:,i),logPPDstore(i), rft_store(:,:,i))
  end do
  
  write(*,*)"Finish computing PPD"

  return 
end subroutine Initialize_pspace

!------------------------------------------------------------
Subroutine AdvanceChain (ichain,T,logPPD)
  implicit none 
  
  Double Precision, intent(out) :: logPPD
  Double Precision, intent(in)  :: T
  Integer, intent(in)          :: ichain
  
  call McMC (ichain,T,logPPD)    ! call users problem specific code for advancing chain
  return
end Subroutine AdvanceChain
!------------------------------------------------------------
!-------------------------------------------------------------------------
!
Subroutine McMC (ichain,T,logPPD)
  use modelspacedata
  use PTmod
  use params
  use mt19937
  implicit none 
  real(8), intent(out) :: logPPD
  real(8), intent(in) :: T
  integer, intent(in) :: ichain
  real(8) :: logQ12, logQ21
  logical :: yn, null_flag
  integer :: tmpk
  real(8) :: tmpv, tmpz, rand, tmpv2
  real(8) :: vp(kmax+2), vs(kmax+2), rho(kmax+2), h(kmax+2)
  integer :: i_type, i_target, i, j, k, iz, isigbin
  integer :: izbin1, izbin2, ivbin, izbin, nlay, ivpvsbin, ihbin, itrc
  real, external :: gasdev

  logPPD = logPPDstore(ichain)
  propmodel(:) = model(:,ichain)
  
  ! count step
  icountMcMC = icountMcMC + 1

  step_count(ichain) = step_count(ichain) + 1

  ! display progress 
  if (ichain == 1 .and. mod(step_count(ichain),100) == 0) then
     write(*,*)step_count(ichain), rank
  end if

  ! make proposal
  null_flag = .false.
  k = int(model(1,ichain))

  ! No need to consider proposal probability as long as 
  ! velocities of newly added layer are extracted from 
  ! the prior probability (see Dosso et al. 2014)
  logQ12 = 0.d0
  logQ21 = 0.d0
  

  i_type = int(grnd() * ntype) + 1
  
  
  if (i_type == itype_birth) then  ! Birth
     tmpk = k + 1
     if (tmpk < kmax) then
        tmpz = zmin + grnd() * del_z
        tmpv = vmin + grnd() * del_v
        tmpv2 = vpvsmin + grnd() * del_vpvs
        propmodel(tmpk+1) = tmpz
        propmodel(tmpk+kmax+1) = tmpv
        propmodel(tmpk+2*kmax+1) = tmpv2
        propmodel(1) = dble(tmpk)
     else
        null_flag = .true.
     end if
  else if (i_type == itype_vs) then  ! velocity perturb
     if (.not. base_flag) then
        i_target = int(grnd()*(k+1)) + 1
        if (i_target == k + 1) then
           i_target = kmax ! bottom half space
        end if
     else 
        i_target = int(grnd()*k) + 1
     end if

     tmpv = propmodel(i_target+kmax+1) + gasdev(iseed) * sdv_v
     if (tmpv >= vmin .and. tmpv <= vmax) then
        propmodel(i_target+kmax+1) = tmpv
     else 
        null_flag = .true.
     end if
  else if (i_type == itype_z) then ! move
     i_target = int(grnd()*k) + 1
     tmpz     = propmodel(1+i_target) + gasdev(iseed) * sdv_z 
     if (tmpz >= zmin .and. tmpz <= zmax) then
        propmodel(1+i_target) = tmpz
     else
        null_flag = .true.
     end if
  else if (i_type == itype_death) then ! Death
     tmpk = k - 1
     if (tmpk >= kmin) then
        i_target = int(grnd()*k) + 1
        do i = i_target, k-1
           propmodel(1+i) = model(2+i,ichain) ! depth
           propmodel(1+kmax+i) = model(2+kmax+i,ichain) ! velocity
           propmodel(1+2*kmax+i) = model(2+2*kmax+i,ichain) ! Vp/Vs
        end do
        propmodel(1+k) = 0.0
        propmodel(1+kmax+k) = 0.0
        propmodel(1+2*kmax+k) = 0.0
        propmodel(1) = dble(tmpk)
     else
        null_flag = .true.
     end if
  else if (i_type == itype_sig) then ! perturb noise sigma
     itrc = grnd() * ntrc + 1
     rand = gasdev(iseed)
     tmpv = propmodel(3*kmax+1+itrc) + rand * sdv_sig
     if (tmpv >= sigmin .and. tmpv <= sigmax) then
        propmodel(3*kmax+1+itrc) = tmpv
     else
        null_flag = .true.
     end if
  else if (i_type == itype_vpvs) then ! perturb Vp/Vs ratio
     if (.not. base_flag) then
        i_target = int(grnd()*(k+1)) + 1
        if (i_target == k + 1) then
           i_target = kmax ! bottom half space
        end if
     else 
        i_target = int(grnd()*k) + 1
     end if
     tmpv = propmodel(i_target+2*kmax+1) + gasdev(iseed) * sdv_vpvs
     if (tmpv >= vpvsmin .and. tmpv <= vpvsmax) then
        propmodel(i_target+2*kmax+1) = tmpv
     else 
        null_flag = .true.
     end if
  else if (i_type == itype_sedvs) then ! sediment Vs
     tmpv = propmodel(3*kmax+ntrc+2) + gasdev(iseed) * sdv_v
     if (tmpv >= sed_vsmin .and. tmpv <= sed_vsmax) then
        propmodel(3*kmax+ntrc+2) = tmpv
     else
        null_flag = .true.
     end if
  else if (i_type == itype_sedvp) then ! sediment Vp
     tmpv = propmodel(3*kmax+ntrc+3) + gasdev(iseed) * sdv_v
     if (tmpv >= sed_vpmin .and. tmpv <= sed_vpmax) then
        propmodel(3*kmax+ntrc+3) = tmpv
     else
        null_flag = .true.
     end if
  else if (i_type == itype_sedh) then ! sediment thickness
     tmpv = propmodel(3*kmax+ntrc+4) + gasdev(iseed) * sdv_z
     if (tmpv >= sed_hmin .and. tmpv <= sed_hmax) then
        propmodel(3*kmax+ntrc+4) = tmpv
     else
        null_flag = .true.
     end if
  else 
     write(*,*)"ERROR"
     stop
  end if
     
  
  
  if (.not. null_flag) then
     call CalclogPPD(int(propmodel(1)), &
          & propmodel(1:npara), logPPD, rft_store(:,:,ichain))
     call PT_McMC_accept(T,logPPDstore(ichain),logQ12,logPPD,logQ21,yn)
     
     if (yn) then
        model(1:npara,ichain) = propmodel(1:npara)
        logPPDstore(ichain) = logPPD
     else
        logPPD = logPPDstore(ichain)
     end if
  end if

  ! record non-tempered chain's sample
  if (T <= 1.d0) then
     n_propose(i_type) = n_propose(i_type) + 1
     if (yn) then
        n_accept(i_type) = n_accept(i_type) + 1
     end if
  end if

  if (icountMcMC > iburnMC .and. T <= 1.d0 .and. &
       & mod(step_count(ichain),nskip) == 0) then
     nmod = nmod + 1

     k = model(1,ichain)
     vs = 0.d0
     h = 0.d0
     
     ! record dimension
     n_dim(k) = n_dim(k) + 1

     ! record sigma
     do itrc = 1, ntrc
        isigbin = int((model(3*kmax+1+itrc,ichain) - sigmin) / del_sigbin) + 1
        n_sig(isigbin,itrc) = n_sig(isigbin,itrc) + 1
     end do

     ! record velocity profile
     call set_struct(int(model(1,ichain)),model(1:npara,ichain),nlay,vp,vs,rho,h)
     do i = 1, k
        izbin = int(model(1+i,ichain) / del_zbin) + 1
        n_int(izbin) = n_int(izbin) + 1
     end do
     
     izbin1 = 1
     do i = 1, nlay
        ivbin = int((vs(i) - vbinmin) / del_vbin) + 1
        if (ivbin < 1) then
           ivbin = 1
        else if (ivbin > n_vbin) then
           ivbin = n_vbin
        end if
        ivpvsbin = int((vp(i)/ vs(i) - vpvsmin) / del_vpvsbin) + 1
        if (ivpvsbin < 1) then
           ivpvsbin = 1
        else if (ivpvsbin > n_vpvsbin) then
           ivpvsbin = n_vpvsbin
        end if
        izbin2 = izbin1 + int(h(i) / del_zbin) 
        if (i == nlay) then
           izbin2 = n_zbin + 1
        end if
        do izbin = izbin1, izbin2 - 1
           if (izbin > n_zbin) cycle
           n_vz(ivbin,izbin) = n_vz(ivbin,izbin) + 1
           n_vpvsz(ivpvsbin,izbin) = n_vpvsz(ivpvsbin,izbin) + 1
        end do
        izbin1 = izbin2
     end do
     
     
     ! record synthetic RFs
     do k = 1, ntrc
        do i = 1, n_smp
           j = int((rft_store(i,k,ichain) - abinmin) / del_abin) + 1
           if (j > n_abin) then
              j = n_abin
           else if (j < 1) then
              j = 1
           end if
           n_amp(i,j,k) = n_amp(i,j,k) + 1
        end do
     end do
  end if
  
  return 
end Subroutine McMC


!-------------------------------------------------------------------------

Subroutine CalclogPPD(k,param,logPPD, rft)
  use params
  use modelspacedata
  use PTmod
  implicit none 
  integer, intent(in) :: k
  real(8), intent(in) :: param(npara)
  real(8), intent(out) :: logPPD
  real(8), intent(out) :: rft(nfft,ntrc)
  real(8) :: misfits(n_smp), phi1(n_smp), phi
  integer :: i, nlay, itrc
  real(8) :: vp(kmax+1), vs(kmax+1), rho(kmax+1), h(kmax+1)
  real(8) :: sig(ntrc), tpre
  real(8) :: llkh
  
  !------------------------------
  ! calculate synthetic RF
  !------------------------------
  rft(:,:) =0.d0
  logPPD = 0.d0
  do itrc = 1, ntrc
     sig(itrc) = param(3*kmax+1+itrc)
  end do
  call set_struct(k,param,nlay,vp,vs,rho,h)
  tpre = -t_win_start

  
  !write(*,*)"========"
  !write(*,*)"test"
  !write(*,*)nfft, delta, 30.0
  !write(*,*)nlay, rayp, 1
  !do i = 1, nlay
  !   write(*,'(4F10.4)')vp(i), vs(i), rho(i), h(i)
  !end do

  call fwd_rf(nlay, nfft, ntrc, rayp, vp, vs, rho, h, rft, & 
       & tpre, a_gus) 
  
  
  llkh = 0.d0
  do itrc = 1, ntrc
     !*** misfits ***
     do i = 1, n_smp
        misfits(i) = rft(i,itrc) - d_obs(i,itrc) 
     end do
     !
     !!*** log-likelihood
     phi1 = matmul(misfits,r_inv(:,:,itrc))
     phi = dot_product(phi1,misfits)

     llkh = llkh - 0.5d0 * phi
     !llkh = llkh - 0.5d0 * phi / (sig(itrc) * sig(itrc)) &
     !     & - dble(n_smp) * log(sig(itrc))
  end do

  !(See Bodin et al. 2012)
  logPPD = -llkh ! must return "negative" log PPD
  
  return 
end Subroutine CalclogPPD

!------------------------------------------------------------

!------------------------------------------------------------
subroutine set_struct(k,param,nlay,vp,vs,rho,h)
  use params
  use PTmod
  implicit none
  integer, intent(in) :: k 
  real(8), intent(in) :: param(npara)
  integer, intent(out) :: nlay
  real(8), intent(out) :: vp(kmax+3), vs(kmax+3), rho(kmax+3), h(kmax+3)
  integer :: i, i_vel_ref
  real(8) :: s1, s2, s3, s4, s5
  real(8) :: tmp_v(kmax), tmp_z(kmax), tmp_vpvs(kmax)
  
  !------------------------------ 
  ! set structure
  !------------------------------ 
  tmp_z(1:k) = param(2:1+k)
  tmp_v(1:k) = param(kmax+2:kmax+1+k)
  tmp_vpvs(1:k) = param(2*kmax+2:2*kmax+1+k)
  if (k > 1) then     
     call hpsort2(k,tmp_z(1:k),tmp_v(1:k),tmp_vpvs(1:k))
  end if
  
  ! *** thickness ***
  nlay = k + 1
  h(1) = tmp_z(1) - zmin
  do i = 2, nlay-1
     h(i) = tmp_z(i) - tmp_z(i-1)
  end do
  h(nlay) = 999.d0

  ! *** physical properties ***
  ! top layer
  i_vel_ref = nint((zmin + 0.5 * h(1)) / dz_vel_ref)+1
  vs(1) = vs_ref(i_vel_ref) + tmp_v(1)
  !vs(1) = tmp_v(1)
  if (vpvs_flag == 0) then
     s1 = vs(1)
     s2 = s1 * s1
     s3 = s2 * s1
     s4 = s3 * s1
     vp(1) = 0.9040 + 2.0947*s1 - 0.8206*s2 + 0.2683*s3 - 0.0251*s4
  else if (vpvs_flag == 1) then
     vp(1) = vp_ref(i_vel_ref)
  else
     vp(1) = vs(1) * tmp_vpvs(1)
  end if
  s1 = vp(1)
  s2 = s1 * s1
  s3 = s2 * s1
  s4 = s3 * s1
  s5 = s4 * s1
  rho(1) = 1.6612*s1 - 0.4721*s2 + 0.0671*s3- 0.0043*s4 &
       & + 0.000106*s5
  do i = 2, nlay-1
     i_vel_ref = nint((tmp_z(i-1) + 0.5 * h(i)) / dz_vel_ref) + 1
     vs(i) = vs_ref(i_vel_ref) + tmp_v(i)
     !vs(i) = tmp_v(i)
     if (vpvs_flag == 0) then
        s1 = vs(i)
        s2 = s1 * s1
        s3 = s2 * s1
        s4 = s3 * s1
        vp(i) = 0.9040 + 2.0947*s1 - 0.8206*s2 + 0.2683*s3 - 0.0251*s4
     else if (vpvs_flag == 1) then
        vp(i) = vp_ref(i_vel_ref)
     else
        vp(i) = vs(i) * tmp_vpvs(i)
     end if
     s1 = vp(i)
     s2 = s1 * s1
     s3 = s2 * s1
     s4 = s3 * s1
     s5 = s4 * s1
     rho(i) = 1.6612*s1 - 0.4721*s2 + 0.0671*s3 - 0.0043*s4 &
          & + 0.000106*s5
  end do
  ! bottom layer
  if (.not. base_flag)  then
     i_vel_ref = nint(0.5d0 * (tmp_z(nlay-1) + zmax) / dz_vel_ref)+1
     vs(nlay) = vs_ref(i_vel_ref) + param(1+2*kmax) 
     
     if (vpvs_flag == 0) then
        s1 = vs(nlay)
        s2 = s1 * s1
        s3 = s2 * s1
        s4 = s3 * s1
        vp(nlay) = 0.9040 + 2.0947*s1 - 0.8206*s2 + 0.2683*s3 - 0.0251*s4
     else if (vpvs_flag == 1) then
        vp(nlay) = vp_ref(i_vel_ref)
     else
        vp(nlay) = vs(nlay) * param(1+3*kmax)
     end if
  else
     vs(nlay) = vs_base
     vp(nlay) = vp_base
  end if
  s1 = vp(nlay)
  s2 = s1 * s1
  s3 = s2 * s1
  s4 = s3 * s1
  s5 = s4 * s1
  rho(nlay) = 1.6612*s1 - 0.4721*s2 + 0.0671*s3 - 0.0043*s4 &
       & + 0.000106*s5

  if (sed_flag == 1) then
     nlay = nlay + 1
     do i = nlay -1, 1, -1
        h(i+1) = h(i)
        rho(i+1) = rho(i)
        vp(i+1) = vp(i)
        vs(i+1) = vs(i)
     end do
     vs(1) = param(3*kmax+ntrc+2)
     vp(1) = param(3*kmax+ntrc+3)
     h(1) = param(3*kmax+ntrc+4)
     rho(1) = 1.8d0

  end if

  if (sdep > 0.d0) then
     nlay = nlay + 1
     do i = nlay -1, 1, -1
        h(i+1) = h(i)
        rho(i+1) = rho(i)
        vp(i+1) = vp(i)
        vs(i+1) = vs(i)
     end do
     vs(1) = -1.d0
     vp(1) = 1.5d0
     rho(1) = 1.d0
     h(1) = sdep
  end if
  
  
  return
end subroutine set_struct
!------------------------------------------------------------

!-------------------------------------------------------------------
!                                               
!       Numerical Recipes random number generator for 
!       a Gaussian distribution
!
! ----------------------------------------------------------------------------



FUNCTION GASDEV(idum)
  
  !  ..Arguments..  integer idum real GASDEV
  
  !     ..Local..
  real v1,v2,r,fac
  real ran3
  
  if (idum.lt.0) iset=0
10 v1=2*ran3(idum)-1
  v2=2*ran3(idum)-1
  r=v1**2+v2**2
  if(r.ge.1.or.r.eq.0) GOTO 10
  fac=sqrt(-2*log(r)/r)
  GASDEV=v2*fac
  
  RETURN
END FUNCTION GASDEV
!
!*****************************************************
!*  Sorts an array RA of length N in ascending order *
!*                by the Heapsort method             *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table RA                   *
!*          RA	  table to be sorted                 *
!*          FA    
!* OUTPUT:                                           *
!*	    RA    table sorted in ascending order    *
!*                                                   *
!* NOTE: The Heapsort method is a N Log2 N routine,  *
!*       and can be used for very large arrays.      *
!*****************************************************         
SUBROUTINE HPSORT(N,RA)
  Double precision RA(N)
  Double precision RRA
  L=N/2+1
  IR=N
  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10 continue
  if(L > 1)then
     L=L-1
     RRA=RA(L)
  else
     RRA=RA(IR)
     RA(IR)=RA(1)
     IR=IR-1
     if(IR.eq.1)then
        RA(1)=RRA
        return
     end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
     if(J < IR)then
        if(RA(J) < RA(J+1))  J=J+1
     end if
     if(RRA < RA(J))then
        RA(I)=RA(J)
        I=J; J=J+J
     else
        J=IR+1
     end if
     goto 20
  end if
  RA(I)=RRA
  goto 10
END SUBROUTINE HPSORT


SUBROUTINE HPSORT2(N,RA, FA, GA)
  Double precision RA(N), FA(N), GA(N)
  Double precision RRA, FFA, GGA
  L=N/2+1
  IR=N
  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10 continue
  if(L > 1)then
     L=L-1
     RRA=RA(L)
     FFA=FA(L)
     GGA=GA(L)
  else
     RRA=RA(IR)
     FFA=FA(IR)
     GGA=GA(IR)
     RA(IR)=RA(1)
     FA(IR)=FA(1)
     GA(IR)=GA(1)
     IR=IR-1
     if(IR.eq.1)then
        RA(1)=RRA
        FA(1)=FFA
        GA(1)=GGA
        return
     end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
     if(J < IR)then
        if(RA(J) < RA(J+1))  J=J+1
     end if
     if(RRA < RA(J))then
        RA(I)=RA(J)
        FA(I)=FA(J)
        GA(I)=GA(J)
        I=J; J=J+J
     else
        J=IR+1
     end if
     goto 20
  end if
  RA(I)=RRA
  FA(I)=FFA
  GA(I)=GGA
  goto 10
END SUBROUTINE HPSORT2
!

!-------------------------------------------------------------------------
!
!     Setuptempladder - User routine to set up temperatures for each chain
!
!     Input:
!        nchains    Integer   : Number of chains
!        Tlow       Double    : Lowest temperature
!        Thigh      Double    : Highest temperature
!
!     Output:
!        Chaintemp  Double array (nchains)    : Temperatures of each chain
!                                               (on this processor)
!           
!     Notes:
!           This utility routine calculates and returns the array Chaintemps.
!           It calculates a random log-uniform set of temperature values
!           between input bounds and put the results in Chaintemps.
!           Chaintemps is written out to file `tlevels'.
!     
!           This routine is provided merely for convenience and the user could set up
!           the arrays Chaintemp in any way desired.
!
!-------------------------------------------------------------------------
!
Subroutine Setuptempladder &
     (nchains,ncool,Tlow,Thigh,Chaintemp)
  use mt19937
  use params, only: iseed
  use PTmod
  
#if defined MPI
  include "mpif.h"
  Integer, dimension(MPI_STATUS_SIZE)           :: status
#endif
  
  Double precision                   :: Chaintemp(nchains)
  Double precision                   :: t,dt,Tlow,Thigh
  Double precision                   :: aval,bval
  Character(len=90)                     filename
  
  ! Selected temperatures randomly using log-uniform distribution  
  
  aval = log(Tlow)
  bval = log(Thigh)
  dx = (bval-aval)/(nchains-1)
  do it=1,nchains
     Chaintemp(it) = exp(aval + grnd()*(bval-aval))
  end do
  

  do it = 1, ncool
     Chaintemp(it) = Tlow    
  end do

  
  allocate(AllTemps(nchains*nproc))
  AllTemps = 0.0
  
#if defined MPI
  ! Send all Temperatures to master for output to file (for diagnostics)
  call MPI_GATHER(Chaintemp,nchains,MPI_DOUBLE_PRECISION,&
       AllTemps,nchains,MPI_DOUBLE_PRECISION,&
       0,MPI_COMM_WORLD,ierror)
  if(rank == 0)then
     filename = trim(dir)//'tlevels'
     open(15,file=filename,status='unknown')
     
     k = nchains*(nproc-1)
     if(k.gt.1)call Hpsort(k,AllTemps(nchains+1)) ! Order the temperatures for neat output
     do it = nchains+1,nproc*nchains
        write(15,*)it-nchains,AllTemps(it)
     end do
     close(15)
  end if
  !deallocate(AllTemps)
#else
  filename = trim(dir)//'tlevels'
  open(15,file=filename,status='unknown')
  do it=1,nchains
     write(15,*)it,Chaintemp(it)
  end do
  AllTemps = Chaintemp
  close(15)
#endif

  return
end Subroutine Setuptempladder
!-------------------------------------------------------------------------
!
!     Setuptempbins - User routine to define temperature bins used for diagnostics
!
!     Input:
!        nbins      Integer   : Number of temperature bins used for diagnostics
!        nchains    Integer   : Number of chains
!        Tlow       Double    : Lowest temperature
!        Thigh      Double    : Highest temperature
!
!     Output:
!        Tbins  Double array (nbins)    : Temperatures of each bin
!        Ttot   Integer array (nbins)   : Number of chains in each bin
!           
!     Notes:
!           Defines a set of temperature bins over which all diagnostic results are 
!           summed. Can be used to look at rates of transition between temperature bins etc.
!
!           Array AllTemps(nchains*nproc) is used to calculate the number of chains
!           across all processors in each temperature bin
!     
!           This routine is provided merely for convenience and the user could set up
!           the arrays Chaintemp in way desired.
!
!-------------------------------------------------------------------------
!
Subroutine Setuptempbins &
     (nbins,Tlow,Thigh,nchains,Tbins,Ttot)
  
  use PTmod
  
  Double precision                   :: Tbins(*)
  Integer                            :: Ttot(*)
  Double precision                   :: dt,Tlow,Thigh
  Double precision                   :: aval,bval
  Character(len=85)                     filename
  
  ! Selected temperatures using log-uniform
  
  aval = log(Tlow)
  bval = log(Thigh)
  dx = (bval-aval)/(nbins-1)
  do it=1,nbins
     Tbins(it) = exp(aval + 0.5*dx + (it-1)*dx)
  end do
  
  if(rank == 0)then
     
     ! Calculate number of temperatures in each temperature bin
     ! for diagnostics output
     do i=1,nbins
        Ttot(i) = 0
     end do
     
     
#if defined MPI
     do it = nchains+1,nproc*nchains
#else
     do it = 1,nchains
#endif
        kbin = nbins
        do i=nbins-1,1,-1
           if(AllTemps(it).le.Tbins(i))then
              kbin = i
           end if
        end do
        Ttot(kbin) = Ttot(kbin) + 1
     end do
     
     ! Write out temperatures on all nodes
     filename = trim(dir)//'Tbins'
     open(15,file=filename,status='unknown')
     write(15,*)'       Index  Temp',&
          '                   No. of chains'
     do it = 1,nbins
        write(15,*)it,Tbins(it),Ttot(it)
     end do
     close(15)
  end if
  
  return
end subroutine Setuptempbins
