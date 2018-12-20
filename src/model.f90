module model
  use params
  use mt19937
  implicit none 
  integer, allocatable :: k(:)
  real(8), allocatable :: vp(:,:), vs(:,:), z(:,:)
  real(8), allocatable :: vp_ref(:), vs_ref(:)
  real(8) :: dz_ref, z_ref_min, z_ref_max
  
contains
  !=====================================================================
  ! generate initial model randomly 
  subroutine init_model()
    implicit none 
    integer :: i, ichain
    
    allocate(z(k_max-1, nchains), vp(k_max, nchains), vs(k_max, nchains))
    allocate(k(nchains))
    
    do ichain = 1, nchains
       k(ichain) = k_min + int(grnd() * (k_max - k_min))
       
       do i = 1, k(ichain)
          z(i, ichain) = z_min + int(grnd() * (z_max - z_min))
       end do
       
       do i = 1, k(ichain) + 1
          vs(i, ichain) = vs_min + grnd() * (vs_max - vs_min)
          vp(i, ichain) = vp_min + grnd() * (vp_max - vp_min)
       end do
    end do
    
    return 
  end subroutine init_model
  !=====================================================================
  subroutine read_ref_model(verb)
    implicit none
    logical, intent(in) :: verb
    integer :: ierr, nref, i
    real(8) :: z_old, z_tmp
    
    ! check if file exists or not
    open(io_ref, file = vel_file, iostat = ierr, status = "old")
    if (ierr /= 0 .and. verb)  then
       write(*,*)"ERROR: cannot open ", trim(vel_file)
       call mpi_finalize(ierr)
       stop
    end if
    
    ! First, obtain # of layers and depth increment
    nref = 0
    dz_ref = -100.d0
    z_old  = -999.d0
    do 
       read(io_ref, *, iostat = ierr)z_tmp
       if (ierr /= 0) exit
       nref = nref + 1
       if (nref >= 3 .and. abs(z_tmp - z_old - dz_ref) > 1.0e-5 &
            & .and. verb) then
          write(*,*)z_tmp, z_old, dz_ref, nref
          write(*,*)"ERROR: Depth increment must be constant in ", &
               & trim(vel_file)
          call mpi_finalize(ierr)
          stop
       end if
       if (nref == 1) then
          z_ref_min = z_tmp
       end if
       dz_ref = z_tmp - z_old
       z_old = z_tmp
    end do
    z_ref_max = z_tmp

    allocate(vp_ref(nref), vs_ref(nref))
    
    ! get reference velocity
    rewind(io_ref)
    do i = 1, nref
       read(io_ref,*) z_tmp, vp_ref(i), vs_ref(i)
    end do
    close(io_ref)
    
    ! Display
    if (verb) then
       write(*,*)
       write(*,*) "--- Reference model ---"
       write(*,*) "# of layer: ", nref
       do i = 1, nref
          if (mod(i,100) == 1) then
             write(*,*)"#", i, z_ref_min + (i - 1) * dz_ref, &
                  & vp_ref(i), vs_ref(i) 
          end if
       end do
    end if


    return
  end subroutine read_ref_model
  
  !=====================================================================
  
  subroutine format_model(ichain, nlay, alpha, beta, h, rho)
    implicit none 
    integer, intent(in) :: ichain
    integer, intent(out) :: nlay
    real(8), intent(out) :: alpha(nlay_max), beta(nlay_max)
    real(8), intent(out) :: h(nlay_max), rho(nlay_max)
    integer :: i, j
    
    i = 0

    ! Ocean
    if (sdep > 0.d0) then
       i = i + 1
       nlay   = k(ichain) + 2
       alpha(i)  = 1.5d0
       beta(i)  = -999.d0
       rho(i) = 1.d0
       h(i)   = sdep
    end if

    ! Finite layer
    do j = 1, k(ichain)
       i = i + 1
       beta(i) = vs(j, ichain)
       if (vp_mode == 1) then
          alpha(i) = vp(j, ichain)
       else
       end if
       rho(i) = vp_to_rho(alpha(i))
       h(i)   = z(j+1, ichain) - z(j, ichain)
    end do

    ! Half space
    i = i + 1
    beta(i) = vs(k(ichain) + 1, ichain)
    if (vp_mode == 1) then
       alpha(i) = vp(j, ichain)
    else

    end if
    rho(i) = vp_to_rho(alpha(i))
    h(i)   = 999.d0

    ! Total number of layers (including ocean and half-space)
    nlay = i
    
    return
  end subroutine format_model
  
  !=====================================================================

  real(8) function vp_to_rho(a1) result(p)
    implicit none 
    real(8), intent(in) :: a1
    real(8) :: a2, a3, a4, a5
    
    a2 = a1 * a1
    a3 = a2 * a1
    a4 = a3 * a1
    a5 = a4 * a1
    
    p = 1.6612 * a1 - 0.4721 * a2 + 0.0671 * a3  - 0.0043 * a4 &
         & + 0.000106 * a5
    ! Empirical relations between density and Vp (Brocher, 2005)


    return 
  end function vp_to_rho

  
end module model
