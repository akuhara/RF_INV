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

module model
  use params
  use mt19937
  use sort
  implicit none 
  integer, allocatable :: k(:)
  real(8), allocatable :: dvp(:,:), dvs(:,:), z(:,:)
  real(8), allocatable :: vp_ref(:), vs_ref(:)
  real(8) :: dz_ref, z_ref_min, z_ref_max
  
contains
  !=====================================================================
  ! generate initial model randomly 
  subroutine init_model(verb)
    use math
    use prior
    implicit none 
    logical, intent(in) :: verb
    integer :: i, ichain
    integer :: nlay
    real(8) :: alpha(nlay_max), beta(nlay_max)
    real(8) :: h(nlay_max), rho(nlay_max)
    real(8), allocatable  :: rft(:,:)
    logical :: is_valid
    
    allocate(z(k_max-1, nchains), dvp(k_max, nchains), dvs(k_max, nchains))
    allocate(k(nchains))
    allocate(rft(nfft, ntrc))

    dvp = 0.d0
    dvs = 0.d0
    z   = 0.d0
    do ichain = 1, nchains
       k(ichain) = k_min + int(grnd() * (k_max - k_min))
       
       do i = 1, k(ichain)
          z(i, ichain) = z_min + grnd() * (z_max - z_min)
       end do
       is_valid = .false.
       do while (.not. is_valid)
          do i = 1, k(ichain)
             !dvs(i, ichain) = dvs_min + grnd() * (dvs_max - dvs_min)
             !dvp(i, ichain) = dvp_min + grnd() * (dvp_max - dvp_min)
             dvs(i, ichain) = gauss() * dvs_prior
             dvp(i, ichain) = gauss() * dvp_prior
          end do
          ! Bottom half-space
          !dvs(k_max, ichain) = dvs_min + grnd() * (dvs_max - dvs_min)
          !dvp(k_max, ichain) = dvp_min + grnd() * (dvp_max - dvp_min)
          dvs(k_max, ichain) = gauss() * dvs_prior
          dvp(k_max, ichain) = gauss() * dvp_prior
          
          call format_model(k(ichain), z(1:k_max-1, ichain), &
               & dvp(1:k_max, ichain), dvs(1:k_max, ichain), &
               & nlay, alpha, beta, rho, h, is_valid)
       end do
       if (verb .and. ichain == 1) then
          write(*,*)
          write(*,*)"--- Initial velocity model ---"
          do i = 1, nlay
             write(*,*)alpha(i), beta(i), rho(i), h(i)
          end do
       end if

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

       ! Reference velocities are adjusted so that they can stay
       ! within the given bounds
       !if (vp_ref(i) + dvp_min < vp_min) then
       !   vp_ref(i) = vp_min - dvp_min
       !end if
       !if (vp_ref(i) + dvp_max > vp_max) then
       !   vp_ref(i) = vp_max - dvp_max
       !end if
       !if (vs_ref(i) + dvs_min < vs_min) then
       !   vs_ref(i) = vs_min - dvs_min
       !end if
       !if (vs_ref(i) + dvs_max > vs_max) then
       !   vs_ref(i) = vs_max - dvs_max
       !end if
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
  
  subroutine format_model(prop_k, prop_z, prop_dvp, prop_dvs, &
       & nlay, alpha, beta, rho, h, is_valid)
    implicit none 
    integer, intent(in) :: prop_k
    real(8), intent(in) :: prop_z(k_max-1), prop_dvp(k_max) 
    real(8), intent(in) :: prop_dvs(k_max)
    integer, intent(out) :: nlay
    real(8), intent(out) :: alpha(nlay_max), beta(nlay_max)
    real(8), intent(out) :: h(nlay_max), rho(nlay_max)
    logical, intent(out) :: is_valid
    real(8) :: zc
    real(8) :: tmp_z(k_max-1), tmp_dvp(k_max), tmp_dvs(k_max)
    integer :: i, j, ki, iz

   
    is_valid = .true.

    tmp_z = prop_z
    tmp_dvs = prop_dvs
    tmp_dvp = prop_dvp

    
    call quick_sort(tmp_z(1:prop_k), 1, prop_k, &
         & tmp_dvp(1:prop_k), tmp_dvs(1:prop_k))

    i = 0
    ! Ocean
    if (sdep > 0.d0) then
       i = i + 1
       alpha(i)  = 1.5d0
       beta(i)  = -999.d0
       rho(i) = 1.d0
       h(i)   = sdep
    end if

    ! Top layer
    i = i + 1
    zc = 0.5d0 * (sdep + tmp_z(1))
    iz = nint((zc - z_ref_min) / dz_ref) + 1
    beta(i) = vs_ref(iz) + tmp_dvs(1)
    if (vp_mode == 1) then
       alpha(i) = vp_ref(iz) + tmp_dvp(1)
    else
       alpha(i) = vp_ref(iz)
    end if
    if (alpha(i) < vp_min .or. alpha(i) > vp_max .or. &
         & beta(i) < vs_min .or. beta(i) > vs_max) then
       is_valid = .false.
    end if
    rho(i) = vp_to_rho(alpha(i))
    h(i)   = tmp_z(1) - sdep

    ! Middle layer
    do j = 2, prop_k
       i = i + 1
       
       zc = 0.5d0 * (tmp_z(j) + tmp_z(j-1))
       iz = nint((zc - z_ref_min) / dz_ref) + 1

       beta(i) = vs_ref(iz) + tmp_dvs(j)
       
       if (vp_mode == 1) then
          alpha(i) = vp_ref(iz) + tmp_dvp(j)
       else
          alpha(i) = vp_ref(iz)
       end if
       if (alpha(i) < vp_min .or. alpha(i) > vp_max .or. &
            & beta(i) < vs_min .or. beta(i) > vs_max) then
          is_valid = .false.
       end if
       rho(i) = vp_to_rho(alpha(i))
       h(i)   = tmp_z(j) - tmp_z(j-1)

    end do

    ! Half space
    i = i + 1
    
    zc = 0.5d0 * (z_max + tmp_z(prop_k))
    iz = nint((zc - z_ref_min) / dz_ref) + 1

    beta(i) = vs_ref(iz) + tmp_dvs(k_max)
    if (vp_mode == 1) then
       alpha(i) = vp_ref(iz) + tmp_dvp(k_max)
    else
       alpha(i) = vp_ref(iz)
    end if
    if (alpha(i) < vp_min .or. alpha(i) > vp_max .or. &
         & beta(i) < vs_min .or. beta(i) > vs_max) then
       is_valid = .false.
    end if
    rho(i) = vp_to_rho(alpha(i))
    h(i)   = 999.d0

    ! Total number of layers (including ocean and half-space)
    nlay = i
    
    
    return 
  end subroutine format_model
  
  !=====================================================================
  subroutine check_model
    
  end subroutine check_model
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
