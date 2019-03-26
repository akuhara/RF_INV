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

module forward
  implicit none 
  !
  real(8), allocatable, public :: flt(:,:)

  !
  complex(kind(0d0)), allocatable, private :: propagator(:, :, :, :, :, :)
  real(8), parameter, private :: pi = 3.1415926535897931d0
  complex(kind(0d0)), parameter, private :: ei = (0.d0, 1.d0)

  
  logical :: is_ray_common

  public init_forward, calc_rf, change_layer_matrix_sol
  private init_filter, calc_seis, &
       & e_inverse, layer_matrix_sol, layer_matrix_liq, &
       & water_level_decon, check_ray, init_propagator
  
contains


  !=====================================================================
  subroutine init_forward(verb)
    implicit none 
    logical, intent(in) :: verb
    
    call init_filter()
    call check_ray(verb)
    call init_propagator()
    
    return 
  end subroutine init_forward
    
  !=====================================================================
  subroutine init_propagator()
    use params
    use model, only: k, dvp, dvs, z, format_model
    implicit none 
    integer :: nlay, ilay, itrc, iomega, nhalf, ichain, ilay0
    real(8) :: omega, df
    real(8) :: alpha(nlay_max), beta(nlay_max), rho(nlay_max), h(nlay_max)
    logical :: is_valid
    complex(kind(0d0)) :: p_mat(4,4)
    
    nhalf = nfft / 2 + 1
    allocate(propagator(4, 4, nlay_max, nhalf, ntrc, nchains))
    propagator = (0.d0, 0.d0)

    df = 1.d0 / (nfft * delta)

    

    do ichain = 1, nchains
       call format_model(k(ichain), z(:,ichain), dvp(:,ichain), &
            & dvs(:,ichain), nlay, alpha, beta, rho, h, is_valid)
       if (beta(1) < 0.d0) then
          ilay0 = 2
       else 
          ilay0 = 1
       end if
       
       do itrc = 1, ntrc
          do iomega = 1, nhalf
             
                
             omega = 2.d0 * pi * df * (iomega - 1)
             if (iomega == 1) then
                omega = 1.0e-5
             end if
             do ilay = ilay0, nlay - 1
                
                call layer_matrix_sol(omega, rho(ilay), alpha(ilay), &
                     & beta(ilay), rayps(itrc), h(ilay), &
                     & propagator(1:4, 1:4, ilay, iomega, itrc, ichain))
                
             end do
          end do
       end do
    end do
    
    return 
  end subroutine init_propagator

  !=====================================================================

  subroutine change_layer_matrix_sol(ilay, ichain, alpha, beta, rho, h)
    use params
    implicit none 
    integer, intent(in) :: ilay, ichain
    real(8), intent(in) :: alpha, beta, rho, h
    real(8) :: df, omega
    integer :: itrc, iomega, nhalf
    
    df = 1.d0 / (nfft * delta)
    nhalf = nfft / 2 + 1
    do itrc = 1, ntrc
       do iomega = 1, nhalf
          omega = 2.d0 * pi * df * (iomega - 1)
          if (iomega == 1) then
             omega = 1.0e-5
          end if
          call layer_matrix_sol(omega, rho, alpha, beta, &
               & rayps(itrc), h, &
               propagator(1:4, 1:4, ilay, iomega, itrc, ichain))
       end do
    end do
    
    return 
  end subroutine change_layer_matrix_sol

  !=====================================================================
  subroutine check_ray(verb)
    use params, only: ntrc, rayps, ipha
    implicit none 
    logical, intent(in) :: verb
    integer :: itrc
    
    is_ray_common = .true.

    if (ntrc > 1) then
       if (verb) then
          write(*,*) "--- check ray parameters ---"
       end if
       do itrc = 2, ntrc
          if (rayps(itrc) /= rayps(1) .or. ipha(itrc) /= ipha(1)) then
             is_ray_common = .false.
          end if
       end do
    end if

    if (verb .and. ntrc > 1) then
       if (is_ray_common) then
          write(*,*) "Ray geometries are common among traces"
          write(*,*) "-> Single FWD mode"
          write(*,*)
       else 
          write(*,*) "Ray geometries are not common among traces"
          write(*,*) "-> Multiple FWD mode"
          write(*,*)
       end if
    end if

    return 
  end subroutine check_ray

  !=====================================================================
    
  subroutine init_filter()
    use params, only: a_gus, ntrc, nfft, delta
    implicit none 
    integer :: itrc, i, nh
    real(8) :: df, fac_norm, omega
    
    nh = nfft / 2 + 1
    df = 1.d0 / (delta * nfft)

    allocate(flt(nh, ntrc))

    do itrc = 1, ntrc
       fac_norm = nfft * a_gus(itrc) * delta / sqrt(pi)
       do i = 1, nh
          omega = (i - 1) * 2.d0 * pi * df
          
          flt(i, itrc) &
               & = exp(-(omega / (2.d0 * a_gus(itrc)))**2) / fac_norm
          
       end do
    end do
    

    return 
  end subroutine init_filter
  
  !=====================================================================
  
  subroutine calc_rf(chain_id, nlay, n, ntrc, rayps, &
       & alpha, beta, rho, h, itype, rft)
    use fftw
    use params, only : delta, a_gus, deconv_mode, t_start, ipha
    implicit none 
    integer, intent(in) :: nlay, n, ntrc, chain_id
    real(8), intent(in) :: rayps(ntrc)
    real(8), intent(in) :: alpha(nlay), beta(nlay), rho(nlay)
    real(8), intent(in) :: h(nlay)
    integer, intent(in) :: itype
    real(8), intent(out) :: rft(n, ntrc)
    complex(kind(0d0)) :: freq_r(n), freq_v(n)
    complex(kind(0d0)) :: rff(n)
    real(8) :: tp, fac_norm
    integer :: nh, itrc, i, npre, j

    nh = n / 2 + 1

    do itrc = 1, ntrc
       
       

       if (itrc == 1 .or. .not. is_ray_common) then
          call calc_seis(itrc, chain_id, itype, nlay, n, rayps(itrc), &
               & ipha(itrc), alpha, beta, rho, h, freq_r, freq_v)
          
          freq_r = conjg(freq_r)
          freq_v = -conjg(freq_v) ! Set upward positive
          
          if (deconv_mode == 1 .and. ipha(itrc) == 1) then
             call water_level_decon(freq_r, freq_v, rff, nh, 0.001d0)
             tp = 0.d0
          else if (deconv_mode == 1 .and. ipha(itrc) == -1) then
             call water_level_decon(freq_v, freq_r, rff, nh, 0.001d0)
             tp = 0.d0
          else
             rff(1:nh) = freq_r(1:nh)
             call direct_P_arrival(nlay, h(1:nlay), alpha(1:nlay), &
                  & rayps(itrc), tp) 
          end if
       end if
       
       ! filter
       cx(1:nh) = rff(1:nh) * flt(1:nh, itrc)
       cx(nh+1:n) = (0.d0, 0.d0)
       
       ! IFFT
       call dfftw_execute(ifft)
       npre = nint((-t_start - tp) / delta)
       
       ! Time shift
       if (ipha(itrc) == 1)  then
          do i = 1, nfft
             j = mod(nfft - npre + i, nfft)
             if (j == 0) then
                j = nfft
             end if
             rft(i, itrc) = rx(j)
          end do
       else
          do i = 1, nfft
             j = mod(nfft + npre - i + 1, nfft)
             if (j == 0) then
                j = nfft
             end if
             rft(i, itrc) = -rx(j)
          end do
       end if
          
       ! Normalizatoin by vertical (only case w/o deconvolution)
       if (deconv_mode == 0) then
          cx(1:nh) = freq_v(1:nh) * flt(1:nh, itrc)
          cx(nh+1:n) = (0.d0, 0.d0)
          call dfftw_execute(ifft)
          fac_norm = maxval(rx)
          rft(:,itrc) = rft(:,itrc) / fac_norm
       end if
       
    end do
    
    return 
  end subroutine calc_rf  
  
  !=====================================================================
  
  subroutine calc_seis(trc_id, chain_id, itype, nlay, npts, &
       & rayp, ipha, alpha, beta, rho, h, ur_freq, uz_freq)
    use params, only: delta, itype_dvs, itype_dvp
    implicit none
    integer, intent(in) :: nlay, ipha, npts, trc_id, chain_id, itype
    real(8), intent(in) :: alpha(nlay), beta(nlay)
    real(8), intent(in) :: rho(nlay), h(nlay)
    real(8), intent(in) :: rayp
    complex(kind(0d0)), intent(out) :: ur_freq(npts), uz_freq(npts)
    real(8) :: omg, domg
    integer :: iomg, nhalf, ilay0, ilay, j, l
    logical :: sea_flag
    complex(kind(0d0)) :: e_inv(4,4), p_prod(4,4), sl(4,4), lq(2,2)
    complex(kind(0d0)) :: denom, a, b, p_mat(4,4)
    
    
    ! Check Whether ocean layer exists
    if (beta(1) < 0) then
       sea_flag = .true.
    else
       sea_flag = .false.
    end if
    nhalf = npts / 2 + 1
    if (sea_flag) then
       ilay0 = 2
    else
       ilay0 = 1
    end if

    domg = 2.d0 * pi /(npts * delta)
    !ur_freq(:) = (0.d0, 0.d0)
    !uz_freq(:) = (0.d0, 0.d0)
    do iomg = 1, nhalf
       omg = dble(iomg - 1) * domg
       if (iomg == 1) then
          omg = 1.0e-5
       end if
       ! Bottom half space
       call e_inverse(omg, rho(nlay), alpha(nlay), beta(nlay), &
            & rayp, e_inv)
       
       ! Clculate production of propagtor matrix of 
       ! all layers
       p_prod = (0.d0, 0.d0)
       do j = 1, 4
          p_prod(j,j) =(1.d0, 0.d0)
       end do
       do ilay = ilay0, nlay - 1

          if (itype == itype_dvp .or. itype == itype_dvs) then
             p_prod = &
                  matmul(propagator(1:4, 1:4, ilay, iomg, trc_id, chain_id), &
                  & p_prod)
          else
             call layer_matrix_sol(omg, rho(ilay), alpha(ilay), &
                  & beta(ilay), rayp, h(ilay), p_mat)
             p_prod = matmul(p_mat, p_prod)
          end if
       end do
       sl = matmul(e_inv, p_prod)

       ! Calculate response using boundary conditions
       if (.not. sea_flag) then
          denom = sl(3,1) * sl(4,2) - sl(3,2) * sl(4,1)
          if (ipha >= 0) then
             ur_freq(iomg) = sl(4,2) / denom
             uz_freq(iomg) = - sl(4,1) / denom
          else 
             ur_freq(iomg) = - sl(3,2) / denom
             uz_freq(iomg) = sl(3,1) / denom
          end if
       else
          call layer_matrix_liq(omg, rho(1), alpha(1), rayp, h(1), lq)
          a = sl(4,2)*lq(1,1) + sl(4,4)*lq(2,1)
          b = sl(3,2)*lq(1,1) + sl(3,4)*lq(2,1)
          if (ipha >= 0) then
             ur_freq(iomg) = a / (a*sl(3,1) - b*sl(4,1))
             uz_freq(iomg) = lq(1,1)*sl(4,1) / (b*sl(4,1)-a*sl(3,1))
          else
             ur_freq(iomg) = - b / (a*sl(3,1) - b*sl(4,1))
             uz_freq(iomg) = - lq(1,1)*sl(3,1) / (b*sl(4,1)-a*sl(3,1))
          end if
       end if
       
    end do
    
    return 
  end subroutine calc_seis

  !=====================================================================
  !------------------------------------------------------------
  ! E_inverse (Aki & Richards, pp. 161, Eq. (5.71))
  !------------------------------------------------------------
  subroutine e_inverse(omega, rho, alpha, beta, p, e_inv)
    implicit none 
    real(8), intent(in) :: omega, p, rho
    real(8), intent(in) :: alpha, beta
    complex(kind(0d0)), intent(out) :: e_inv(4,4)
    real(8) :: eta, xi, bp
    
    e_inv(:,:) = (0.d0, 0.d0)
    eta = sqrt(1.d0/(beta*beta) - p*p)
    xi  = sqrt(1.d0/(alpha*alpha) - p*p)
    bp = 1.d0 - 2.d0*beta*beta*p*p
    
    e_inv(1,1) = beta*beta*p/alpha
    e_inv(1,2) = bp/(2.d0*alpha*xi)
    e_inv(1,3) = -p/(2.d0*omega*rho*alpha*xi) * ei
    e_inv(1,4) = -1.d0/(2.d0*omega*rho*alpha) * ei
    e_inv(2,1) = bp / (2.d0*beta*eta)
    e_inv(2,2) = -beta*p
    e_inv(2,3) = -1.0/(2.d0*omega*rho*beta) * ei
    e_inv(2,4) = p/(2.d0*omega*rho*beta*eta) * ei
    e_inv(3,1) = e_inv(1,1)
    e_inv(3,2) = - e_inv(1,2)
    e_inv(3,3) = - e_inv(1,3)
    e_inv(3,4) = e_inv(1,4)
    e_inv(4,1) = e_inv(2,1)
    e_inv(4,2) = - e_inv(2,2)
    e_inv(4,3) = - e_inv(2,3)
    e_inv(4,4) = e_inv(2,4)
    
    return 
  end subroutine e_inverse
  
  !------------------------------------------------------------
  ! layer_matrix (Aki & Richards, pp. 398, Eq. (3) in Box 9.1)
  !------------------------------------------------------------
  subroutine layer_matrix_sol(omega, rho, alpha, beta, p, z, p_mat)
    implicit none
    real(8), intent(in) :: alpha, beta
    real(8), intent(in) :: omega, rho, p, z
    complex(kind(0d0)), intent(out) :: p_mat(4,4)
    real(8) :: eta,xi,beta2,p2,bp,cos_xi,cos_eta,sin_xi,sin_eta
    
    beta2 = beta*beta
    p2 =p*p
    bp = 1.d0 -2.d0*beta2*p2
    eta = sqrt(1.d0/(beta2) - p2)
    xi  = sqrt(1.d0/(alpha*alpha) - p2)
    cos_xi = cos(omega*xi*z)
    cos_eta = cos(omega*eta*z)
    sin_xi = sin(omega*xi*z)
    sin_eta = sin(omega*eta*z)
    

    p_mat(1,1) = 2.d0*beta2*p2*cos_xi + bp*cos_eta
    p_mat(2,1) = p*( 2.d0*beta2*xi*sin_xi - bp/eta*sin_eta ) * ei
    p_mat(3,1) = omega*rho*( -4.d0*beta2*beta2*p2*xi*sin_xi - bp*bp/eta*sin_eta )
    p_mat(4,1) = 2.d0*omega*beta2*rho*p*bp*( cos_xi - cos_eta ) * ei
    p_mat(1,2) = p*( bp/xi*sin_xi - 2.d0*beta2*eta*sin_eta ) * ei
    p_mat(2,2) = bp*cos_xi + 2.d0*beta2*p2*cos_eta
    p_mat(3,2) = p_mat(4,1)
    p_mat(4,2) = -omega*rho*( bp*bp/xi*sin_xi + 4.d0*beta2*beta2*p2*eta*sin_eta  )    
    p_mat(1,3) = (p2/xi*sin_xi + eta*sin_eta)/(omega*rho)
    p_mat(2,3) = p*(-cos_xi + cos_eta)/(omega*rho) * ei  
    p_mat(3,3) = p_mat(1,1)
    p_mat(4,3) = p_mat(1,2)  
    p_mat(1,4) = p_mat(2,3)
    p_mat(2,4) = (xi*sin_xi + p2/eta*sin_eta)/(omega*rho)
    p_mat(3,4) = p_mat(2,1)
    p_mat(4,4) = p_mat(2,2)
    
    return 
  end subroutine layer_matrix_sol
  
  !------------------------------------------------------------
  subroutine layer_matrix_liq(omega, rho, alpha, p, z, p_mat)
    implicit none 
    real(8), intent(in) :: omega, rho, alpha, p, z
    complex(kind(0d0)), intent(out) :: p_mat(2,2)
    real(8) :: xi, cos_xi, sin_xi, g
    
    
    xi  = sqrt(1.d0/(alpha*alpha) - p * p)
    cos_xi = cos(omega*xi*z)
    sin_xi = sin(omega*xi*z)
    g = rho * omega / xi
    
    p_mat(1,1) = cos_xi
    p_mat(1,2) = sin_xi / g
    p_mat(2,1) = -g * sin_xi
    p_mat(2,2) = cos_xi
    
    return 
  end subroutine layer_matrix_liq
  
  !=====================================================================
  
  !-----------------------------------------------------------------------
  subroutine water_level_decon(y,x,z,n,pcnt) ! z = y / x
    implicit none 
    integer, intent(in) :: n
    complex(kind(0d0)), intent(in) :: x(n), y(n)
    complex(kind(0d0)), intent(out) :: z(n)
    real(8), intent(in) :: pcnt
    integer :: i 
    real(8) :: wlvl
    real(8) :: amp(n)
    
    do i = 1, n
       amp(i) = x(i) * conjg(x(i))
    end do
    wlvl = pcnt * maxval(amp)
    
    
    do i = 1, n
       z(i) = y(i) * conjg(x(i))/ max(amp(i), wlvl)
    end do

    !z(1:n) = y(1:n) / x(1:n)
    
    return
  end subroutine water_level_decon
  
  !---------------------------------------------------------------------

  subroutine direct_P_arrival(nlay, h, v, rayp, t)
    implicit none 
    integer, intent(in) :: nlay
    real(8), intent(in) :: rayp, h(nlay), v(nlay)
    real(8), intent(out) :: t
    integer :: i, i0
    
    t = 0.d0
    if (v(nlay) >= 0.d0) then
       i0 = 1
    else
       i0 = 2 ! sea water
    end if
    do i = 1, nlay-1
       t = t + h(i) * sqrt(1.d0 / (v(i)* v(i)) - rayp * rayp)
    end do

    return 
  end subroutine direct_P_arrival
    
  
end module forward
