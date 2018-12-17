module params
  implicit none 

  ! constant
  integer, parameter :: clen_max = 200
  integer, parameter :: io_param = 10
  
  ! iteration
  integer :: nburn, niter, ncorr
  
  ! parallel tempering
  integer :: nchains, ncool
  real(8) :: t_high
  
  ! random number seed
  integer :: iseed
  
  ! observation
  integer :: ntrc
  character(clen_max), allocatable :: obs_files(:), &
       & err_files(:), acf_files(:)
  real(8), allocatable :: rayps(:), a_gus(:)
  real(8) :: delta, t_start, t_end
  real(8) :: sdep
  
  ! reference velocity
  character(clen_max) :: vel_file
  
  ! inversion setting
  integer :: vp_mode
  
  ! prior
  integer :: k_min, k_max
  real(8) :: z_min, z_max
  real(8) :: vs_min, vs_max, vpvs_min, vpvs_max
  real(8) :: dev_vs_prior, dev_vpvs_prior
  
  ! proposal
  real(8) :: dev_z, dev_vs, dev_vpvs

  ! output
  integer :: nbin_z, nbin_vs, nbin_vpvs, nbin_amp
  real(8) :: amp_min, amp_max

  !=====================================================================
  
contains
  subroutine get_params(verb, param_file)
    logical, intent(in) :: verb
    character(*), intent(in) :: param_file
    character(100) :: line
    integer :: ierr
    
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

    if (verb) then
       write(*,*)"Nburn = ", nburn
       write(*,*)"Niter = ", niter
       write(*,*)"Ncorr = ", ncorr
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
  
end module params
