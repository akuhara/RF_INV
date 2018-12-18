module params
  implicit none 

  ! constant
  integer, parameter :: clen_max = 200
  integer, parameter :: io_param = 10
  integer, parameter :: npts_max = 2000
  
  ! iteration
  integer :: nburn, niter, ncorr
  
  ! parallel tempering
  integer :: nchains, ncool
  real(8) :: t_high
  
  ! random number seed
  integer :: iseed
  
  ! observation
  integer :: ntrc, nsmp
  character(clen_max), allocatable :: obs_files(:)
  real(8), allocatable :: rayps(:), a_gus(:)
  real(8), allocatable :: obs(:,:)
  real(8) :: delta, t_start, t_end
  real(8) :: sdep

  
  ! reference velocity
  character(clen_max) :: vel_file
  
  ! inversion setting
  integer :: vp_mode
  
  ! prior
  integer :: k_min, k_max
  real(8) :: z_min, z_max
  real(8) :: vs_min, vs_max, vp_min, vp_max
  real(8) :: dev_vs_prior, dev_vp_prior
  
  ! proposal
  real(8) :: dev_z, dev_vs, dev_vp

  ! output
  integer :: nbin_z, nbin_vs, nbin_vp, nbin_amp
  real(8) :: amp_min, amp_max

  !=====================================================================
  
contains
  subroutine get_params(verb, param_file)
    logical, intent(in) :: verb
    character(*), intent(in) :: param_file
    character(100) :: line
    integer :: ierr, itrc
    
    open(io_param, file = param_file, status = "old", iostat = ierr)
    if (ierr /= 0) then
       write(0,*)"ERROR: cannot open : params.in"
       stop
    end if
    
    call get_line(io_param, line)
    read(line,*) nburn

    call get_line(io_param, line)
    read(line,*) niter

    call get_line(io_param, line)
    read(line,*) ncorr
    
    call get_line(io_param, line)
    read(line,*) nchains

    call get_line(io_param, line)
    read(line,*) ncool

    call get_line(io_param, line)
    read(line,*) t_high

    call get_line(io_param, line)
    read(line,*) iseed

    call get_line(io_param, line)
    read(line,*) ntrc
    
    call allocate_trace_num()
    
    do itrc = 1, ntrc
       call get_line(io_param, line)
       read(line,*) rayps(itrc)
    end do
    
    do itrc = 1, ntrc
       call get_line(io_param, line)
       read(line,*) a_gus(itrc)
    end do
    
    do itrc = 1, ntrc
       call get_line(io_param, line)
       read(line,*) obs_files(itrc)
    end do

    call get_line(io_param, line)
    read(line,*) t_start, t_end

    call get_line(io_param, line)
    read(line,*) sdep

    call get_line(io_param, line)
    read(line,*) vel_file

    call get_line(io_param, line)
    read(line,*) vp_mode
    
    call get_line(io_param, line)
    read(line,*) k_min, k_max

    call get_line(io_param, line)
    read(line,*) z_min, z_max

    call get_line(io_param, line)
    read(line,*) vs_min, vs_max

    call get_line(io_param, line)
    read(line,*) vp_min, vp_max

    call get_line(io_param, line)
    read(line,*) dev_vs_prior

    call get_line(io_param, line)
    read(line,*) dev_vp_prior

    call get_line(io_param, line)
    read(line,*) dev_z

    call get_line(io_param, line)
    read(line,*) dev_vs
    
    call get_line(io_param, line)
    read(line,*) dev_vp

    call get_line(io_param, line)
    read(line,*) nbin_z
    
    call get_line(io_param, line)
    read(line,*) nbin_vs
    
    call get_line(io_param, line)
    read(line,*) nbin_vp

    call get_line(io_param, line)
    read(line,*) nbin_amp

    call get_line(io_param, line)
    read(line,*) amp_min, amp_max

    if (verb) then
       write(*,*)"--- Parameters --- "
       write(*,*)"# of iteration in burn-in              : ", nburn
       write(*,*)"# of iteration after burn-in           : ", niter
       write(*,*)"# of iteration per sample              : ", ncorr
       write(*,*)"# of chains per processor              : ", nchains
       write(*,*)"# of non-tempered chains per processor : ", ncool
       write(*,*)"Highest temperature                    : ", t_high
       write(*,*)"Random seed number                     : ", iseed
       write(*,*)"# of observed RFs                      : ", ntrc
       do itrc = 1, ntrc
          write(*,*)"Ray parameter                          : ", &
               & rayps(itrc), "s/km"
       end do
       do itrc = 1, ntrc
          write(*,*)"Gaussian parameter                     : ", &
               & a_gus(itrc)
       end do
       do itrc = 1, ntrc
          write(*,*)"Observed RF file                       : ", &
               & trim(obs_files(itrc))
       end do
       write(*,*)"Start/End time                         : ", t_start, &
            & t_end
       write(*,*)"Station depth                          : ", sdep
       write(*,*)"Reference velocity file                : ", trim(vel_file)
       write(*,*)"Vp mode (0: Fixed, 1: Solved)          : ", vp_mode
       write(*,*)"Min./Max. # of interfaces              : ", k_min, k_max
       write(*,*)"Min./Max. of interface depth           : ", z_min, z_max
       write(*,*)"Min./Max. of Vs prior                  : ", vs_min, vs_max
       write(*,*)"Min./Max. of Vp prior                  : ", vp_min, vp_max
       write(*,*)"Standard deviation for Vs prior        : ", dev_vs_prior
       write(*,*)"Standard deviation for Vp prior        : ", dev_vp_prior
       write(*,*)"Standard deviation for depth proposal  : ", dev_z
       write(*,*)"Standard deviation for Vs proposal     : ", dev_vs
       write(*,*)"Standard deviation for Vp proposal     : ", dev_vp
       write(*,*)"# of bins for depth                    : ", nbin_z
       write(*,*)"# of bins for Vs                       : ", nbin_vs
       write(*,*)"# of bins for Vp                       : ", nbin_vp
       write(*,*)"# of bins for amplitudes               : ", nbin_amp
       write(*,*)"Min./Max. amplitudes to be displayed   : ", amp_min, amp_max
       
    end if
    

    close(io_param)
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
    
    allocate(rayps(ntrc), a_gus(ntrc))
    allocate(obs_files(ntrc), obs(npts_max, ntrc))

    
    return 
  end subroutine allocate_trace_num
  
end module params
