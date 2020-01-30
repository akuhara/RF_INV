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

module params
  implicit none 

  ! Variables shown here are all for pulic use
  
  ! constant
  integer, parameter :: clen_max = 200
  integer, parameter :: io_param = 10, io_obs = 20, io_ref = 30
  integer, parameter :: io_nk = 40, io_syn = 50, io_vpz = 60
  integer, parameter :: io_vsz = 70, io_vpvsz = 80, io_z = 90
  integer, parameter :: io_sig = 100
  integer, parameter :: io_copy = 110, io_lkhd = 120
  integer, parameter :: io_vp_mean = 130, io_vs_mean = 140
  integer, parameter :: io_vpvs_mean = 150
  

  integer, parameter :: npts_max = 2000, nlay_max = 200

  ! Output directory
  character(clen_max) :: out_dir

  ! iteration
  integer :: nburn, niter, ncorr
  
  ! parallel tempering
  integer :: nchains, ncool
  real(8) :: t_high
  real(8), allocatable :: temps(:)
  
  ! random number seed
  integer :: iseed
  
  ! observation
  integer :: ntrc, nsmp
  integer, allocatable :: ipha(:)
  character(clen_max), allocatable :: obs_files(:)
  real(8), allocatable :: rayps(:)
  real(8), allocatable :: obs(:,:)
  real(8) :: delta, t_start, t_end
  real(8) :: sdep, bdep

  ! Receiver function
  integer :: nfft, deconv_mode
  real(8), allocatable :: a_gus(:)


  ! reference velocity
  character(clen_max) :: vel_file
  
  ! inversion setting
  integer :: vp_mode
  integer, allocatable :: sig_mode(:)
  
  
  ! prior
  integer :: k_min, k_max
  real(8) :: z_min, z_max, h_min

  integer :: prior_mode
  real(8) :: dvs_prior, dvp_prior
  real(8), allocatable :: sig_min(:), sig_max(:)
  
  ! proposal
  real(8) :: dev_z, dev_dvs, dev_dvp, dev_sig

  ! output
  integer :: nbin_z, nbin_vs, nbin_vp, nbin_amp, nbin_sig, nbin_vpvs
  real(8) :: amp_min, amp_max, vp_min, vp_max, vs_min, vs_max
  real(8) :: vpvs_min, vpvs_max

  !=====================================================================
  
contains
  subroutine get_params(verb, param_file)
    logical, intent(in) :: verb
    character(*), intent(in) :: param_file
    character(clen_max) :: line, out_file
    integer :: ierr, itrc, ierr2
    
    open(io_param, file = param_file, status = "old", iostat = ierr)
    if (ierr /= 0) then
       write(0,*)"ERROR: cannot open : params.in"
       stop
    end if

    call get_line(io_param, line)
    read(line,*) out_dir
    
    out_file = trim(out_dir) // "/" // "params.in.copy"
    open(io_copy, file = out_file, status = "unknown", iostat = ierr)
    if (ierr /= 0) then
       write(0,*)"ERROR: cannot create ", trim(out_file)
       call mpi_finalize(ierr)
       stop
    end if

    write(io_copy, *) trim(out_dir)

    call get_line(io_param, line)
    read(line,*) nburn
    write(io_copy, *) nburn

    call get_line(io_param, line)
    read(line,*) niter
    write(io_copy, *) niter
    
    call get_line(io_param, line)
    read(line,*) ncorr
    write(io_copy, *) ncorr
    
    call get_line(io_param, line)
    read(line,*) nchains
    write(io_copy, *) nchains
    
    call get_line(io_param, line)
    read(line,*) ncool
    write(io_copy, *) ncool
    
    call get_line(io_param, line)
    read(line,*) t_high
    write(io_copy, *) t_high
    
    call get_line(io_param, line)
    read(line,*) iseed
    write(io_copy, *) iseed
    
    call get_line(io_param, line)
    read(line,*) ntrc
    write(io_copy, *) ntrc
    
    call allocate_trace_num()
    
    do itrc = 1, ntrc
       call get_line(io_param, line)
       read(line,*) rayps(itrc)
       write(io_copy, *) rayps(itrc)
    end do
    
    do itrc = 1, ntrc
       call get_line(io_param, line)
       read(line,*) a_gus(itrc)
       write(io_copy, *) a_gus(itrc)
    end do
    
    do itrc = 1, ntrc
       call get_line(io_param, line)
       read(line,*) ipha(itrc)
       write(io_copy, *) ipha(itrc)
    end do

    call get_line(io_param, line)
    read(line,*) nfft
    write(io_copy, *) nfft
    
    do itrc = 1, ntrc
       call get_line(io_param, line)
       read(line,*) obs_files(itrc)
       write(io_copy, *) trim(obs_files(itrc))
    end do

    call get_line(io_param, line)
    read(line,*) t_start, t_end
    write(io_copy, *) t_start, t_end
    
    call get_line(io_param, line)
    read(line,*)deconv_mode
    write(io_copy, *) deconv_mode
    if (deconv_mode /= 0 .and. deconv_mode /= 1) then
       write(0,*) "ERROR: deconv_mode must be either 0 or 1"
       call mpi_finalize(ierr)
       stop
    end if
    
    
    call get_line(io_param, line)
    read(line,*,iostat=ierr) sdep
    bdep = 0.d0
    if (ierr /=0) then
       read(line,*,iostat=ierr2)sdep, bdep
       if (ierr2 /= 0) then
          write(0,*)"ERROR: while reading SEA_DEP" // &
               & " (BORE_HOLE_DEP)"
          call mpi_finalize(ierr)
          stop
       end if
    end if
       
    write(io_copy, *) sdep

    call get_line(io_param, line)
    read(line,*) vel_file
    write(io_copy, *) trim(vel_file)

    call get_line(io_param, line)
    read(line,*) vp_mode
    write(io_copy, *) vp_mode

    call get_line(io_param, line)
    read(line,*) k_min, k_max
    write(io_copy, *) k_min, k_max

    call get_line(io_param, line)
    read(line, *) z_min, z_max
    write(io_copy, *) z_min, z_max
    
    call get_line(io_param, line)
    read(line, *) h_min
    write(io_copy, *) h_min

    call get_line(io_param, line)
    read(line, *) prior_mode
    write(io_copy, *) prior_mode

    call get_line(io_param, line)
    read(line,*) dvs_prior
    write(io_copy, *) dvs_prior

    call get_line(io_param, line)
    read(line,*) dvp_prior
    write(io_copy, *) dvp_prior

    do itrc = 1, ntrc
       call get_line(io_param, line)
       read(line,*) sig_min(itrc), sig_max(itrc)
       write(io_copy, *) sig_min(itrc), sig_max(itrc)
       if (sig_max(itrc) - sig_min(itrc) > 1.0e-5) then
          sig_mode(itrc) = 1
          if (verb) then
             write(*,*)"Sigma is solved for trace", itrc
          end if
       else
          sig_mode(itrc) = 0
          if (verb) then
             write(*,*) "Sigma is fixed at ", sig_min(itrc) &
                  & , " for trace ", itrc
          end if
       end if
    end do
    
    call get_line(io_param, line)
    read(line,*) dev_z
    write(io_copy, *) dev_z

    call get_line(io_param, line)
    read(line,*) dev_dvs
    write(io_copy, *) dev_dvs
    
    call get_line(io_param, line)
    read(line,*) dev_dvp
    write(io_copy, *) dev_dvp

    call get_line(io_param, line)
    read(line,*) dev_sig
    write(io_copy, *) dev_sig

    call get_line(io_param, line)
    read(line,*) nbin_z
    write(io_copy, *) nbin_z
    
    call get_line(io_param, line)
    read(line,*) nbin_vs
    write(io_copy, *) nbin_vs
    
    call get_line(io_param, line)
    read(line,*) nbin_vp
    write(io_copy, *) nbin_vp

    call get_line(io_param, line)
    read(line, *) nbin_vpvs
    write(io_copy, *)nbin_vpvs

    call get_line(io_param, line)
    read(line,*) nbin_sig
    write(io_copy, *) nbin_sig

    call get_line(io_param, line)
    read(line,*) nbin_amp
    write(io_copy, *) nbin_amp

    call get_line(io_param, line)
    read(line,*) amp_min, amp_max
    write(io_copy, *) amp_min, amp_max

    call get_line(io_param, line)
    read(line,*) vp_min, vp_max
    write(io_copy, *) vp_min, vp_max

    call get_line(io_param, line)
    read(line,*) vs_min, vs_max
    write(io_copy, *) vs_min, vs_max
    
    call get_line(io_param, line)
    read(line, *) vpvs_min, vpvs_max
    write(io_copy, *) vpvs_min, vpvs_max

    if (verb) then
       write(*,*)"--- Parameters --- "
       write(*,*)"OUT_DIR: ", trim(out_dir)
       write(*,*)"N_BURN: ", nburn
       write(*,*)"N_ITER: ", niter
       write(*,*)"N_CORR: ", ncorr
       write(*,*)"N_CHAINS: ", nchains
       write(*,*)"N_COOL: ", ncool
       write(*,*)"T_HIGH: ", t_high
       write(*,*)"I_SEED: ", iseed
       write(*,*)"N_TRC: ", ntrc
       do itrc = 1, ntrc
          write(*,*)"RAYP: ", rayps(itrc), "s/km"
       end do
       do itrc = 1, ntrc
          write(*,*)"A_GAUSS: ", a_gus(itrc)
       end do
       do itrc = 1, ntrc
          write(*,*)"I_PHA: ", ipha(itrc)
       end do
       write(*,*)"N_FFT: ", nfft
       do itrc = 1, ntrc
          write(*,*)"OBS_FILES: ", trim(obs_files(itrc))
       end do
       write(*,*)"T_START, T_END: ", t_start, t_end
       write(*,*)"DEONV_MODE: ", deconv_mode
       write(*,*)"SEA_DEP: ", sdep
       write(*,*)"VEL_FILE: ", trim(vel_file)
       write(*,*)"VP_MODE: ", vp_mode
       write(*,*)"K_MIN K_MAX: ", k_min, k_max
       write(*,*)"Z_MIN Z_MAX: ", z_min, z_max
       write(*,*)"H_MIN: ", h_min        
       write(*,*)"PRIOR_TYPE: ", prior_mode
       write(*,*)"DEV_DVS_PRIOR: ", dvs_prior
       write(*,*)"DEV_DVP_PRIOR: ", dvp_prior
       do itrc = 1, ntrc
          write(*,*)"SIG_MIN SIG_MAX: ", sig_min(itrc), sig_max(itrc)
       end do
       write(*,*)"STEP_SIZE_Z: ", dev_z
       write(*,*)"STEP_SIZE_DVS: ", dev_dvs
       write(*,*)"STEP_SIZE_DVP: ", dev_dvp
       write(*,*)"STEP_SIZE_SIG: ", dev_sig
       write(*,*)"N_BIN_Z: ", nbin_z
       write(*,*)"N_BIN_VS: ", nbin_vs
       write(*,*)"N_BIN_VP: ", nbin_vp
       write(*,*)"N_VIN_VPVS: ", nbin_vpvs
       write(*,*)"N_BIN_SIG: ", nbin_sig
       write(*,*)"N_BIN_AMP: ", nbin_amp
       write(*,*)"AMP_MIN AMP_MAX: ", amp_min, amp_max
       write(*,*)"VP_MIN VP_MAX: ", vp_min, vp_max
       write(*,*)"VS_MIN VS_MAX: ", vs_min, vs_max
       write(*,*)"VPVS_MIN VPVS_MAX: ", vpvs_min, vpvs_max
    end if
    close(io_param)
    close(io_copy)
    return 
  end subroutine get_params
  
  !=====================================================================
  
  subroutine get_line(i_unit,line)
    integer, intent(in) :: i_unit
    character(100), intent(out) :: line
    do 
       read(i_unit,'(a)')line
       line = adjustl(line)
       if (line(1:1) == "#") then
          cycle
       else
          return
       end if
    end do
    return 
  end subroutine get_line
  
  !=====================================================================
  
  subroutine allocate_trace_num()
    implicit none 

    allocate(rayps(ntrc), a_gus(ntrc), ipha(ntrc))
    allocate(obs_files(ntrc), obs(npts_max, ntrc))
    allocate(sig_min(ntrc), sig_max(ntrc))
    allocate(sig_mode(ntrc))
    
    return 
  end subroutine allocate_trace_num
  
  !=====================================================================
  
  subroutine read_obs(verb)
    implicit none 
    logical, intent(in) :: verb
    integer :: itrc, ierr, npts, it1, it2, it
    real :: delta4, t_beg4, tmp4(npts_max)
    character(clen_max) :: ofile
    
    
    if (verb) then
       write(*,*)
       write(*,*)"--- Reading observed data ---"
    end if
    
    
    do itrc = 1, ntrc
       open(io_obs, file = obs_files(itrc), status = "old", &
            & access = "direct", recl = 4, iostat = ierr)
       if (ierr /= 0 ) then
          if (verb) write(*,*) "ERROR: cannot open ", &
               & trim(obs_files(itrc))
          call mpi_finalize(ierr)
          stop
       end if
       
       read(io_obs, rec = 1) delta4
       read(io_obs, rec = 6) t_beg4
       read(io_obs, rec = 80) npts
       it1 = nint((t_start - t_beg4) / delta4) + 1
       it2 = nint((t_end - t_beg4) / delta4) + 1
       nsmp = it2 - it1 + 1
       delta = dble(delta4)
       
       
       do it = 1, nsmp
          read(io_obs, rec = 158 + it + it1 - 1) tmp4(it)
       end do
       obs(1:nsmp, itrc) = tmp4(1:nsmp)
       close(io_obs)
       
       write(ofile, '(A5,I2.2)')"input", itrc
       
       open(io_obs, file = ofile, status = "unknown")
       do it = 1, nsmp
          write(io_obs,*)(it - 1) * delta + t_start, obs(it, itrc)
       end do
       close(io_obs)
       
       if (verb) then
          write(*,*) "Finish reading ", trim(obs_files(itrc))
          write(*,*) "  output -> ", trim(ofile)
       end if
    end do
    
    return 
  end subroutine read_obs

end module params
