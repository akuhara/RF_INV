module forward
  implicit none 
  
  real(8), allocatable, private :: flt(:,:)
  real(8), parameter, private :: pi = 3.1415926535897931
  complex(kind(0d0)), parameter, private :: ei = (0.d0, 1.d0)
  
contains
  
  !=====================================================================

  subroutine fwd_rf(nlay, n, ntrc, rayps, alpha, beta, rho, h, rft)
    use fftw
    implicit none 
    
    ! I/O variables
    integer, intent(in) :: nlay, n, ntrc
    real(8), intent(in) :: rayps(ntrc)
    real(8), intent(in) :: alpha(nlay), beta(nlay), rho(nlay)
    real(8), intent(in) :: h(nlay)
    real(8), intent(out) :: rft(n, ntrc)
    !real(8) :: dtmp(n)
    complex(kind(0d0)) :: freq_r(n), freq_v(n)
    complex(kind(0d0)) :: rff(n)
    integer :: nh, itrc, i
    

    nh = n / 2 + 1
    
    do itrc = 1, ntrc
       call fwd_seis(nlay, n, rayps(itrc), 1, alpha, beta, rho, h, &
            & freq_r, freq_v)
       freq_r = conjg(freq_r)
       freq_v = -conjg(freq_v)
       
       call water_level_decon(freq_r, freq_v , rff, nh, 0.005d0)

       ! filter
       rff(1:nh) = rff(1:nh) * flt(1:nh, itrc)
       
       ! Radial
       cx(1:nh) = rff(1:nh)
       cx(nh+1:n) = 0.d0
       call dfftw_execute(ifft)
       rft(1:n, itrc) = rx(1:n)
       
    end do
    

    return 
  end subroutine fwd_rf
  
  !=====================================================================
  
  subroutine init_filter()
    use params, only: a_gus, ntrc, nfft, delta
    implicit none 
    integer :: itrc, i, nh
    real(8) :: df, fac_norm

    nh = nfft / 2 + 1
    df = 1.d0 / (delta * nfft)

    allocate(flt(nh, ntrc))

    
    do itrc = 1, ntrc
       fac_norm = nfft * a_gus(itrc) * delta / sqrt(pi)
       do i = 1, nh
          flt(i, itrc) &
               & = exp(-(2.d0 * pi * (i - 1) * df &
               & / (2.d0 * a_gus(itrc)))**2) / &
               & fac_norm
       end do
    end do
    

    return 
  end subroutine init_filter
  
  !=====================================================================
  
  subroutine fwd_seis (nlay, npts, rayp, ipha, &
       & alpha, beta, rho, h, ur_freq, uz_freq)
    use params, only: delta
    implicit none
    integer, intent(in) :: nlay, ipha, npts
    real(8), intent(in) :: alpha(nlay), beta(nlay)
    real(8), intent(in) :: rho(nlay), h(nlay)
    real(8), intent(in) :: rayp
    complex(kind(0d0)), intent(out) :: ur_freq(npts), uz_freq(npts)
    real(8) :: omg, domg
    integer :: i, nhalf
    logical :: sea_flag
    
    !do i = 1, nlay
    !   write(*,*)alpha(i), beta(i), rho(i), h(i)
    !end do

    if (beta(1) < 0) then
       sea_flag = .true.
    else
       sea_flag = .false.
    end if
    nhalf = npts / 2 + 1
    
    domg = 2.d0 * pi /(npts * delta)
    ur_freq(1) = (0.d0, 0.d0)
    uz_freq(1) = (0.d0, 0.d0)
    do i = 2, nhalf
       omg = dble(i-1) * domg
       call haskell(omg, ur_freq(i), uz_freq(i), rayp, nlay, ipha, &
            & alpha, beta, rho, h, sea_flag)
    end do
    
    return 
  end subroutine fwd_seis
  
  !=====================================================================

  subroutine haskell(omega, ur, uz, p, nl, ipha, alpha, beta, rho, h, &
       & sea_flag)
    implicit none 
    
    ! I/O variables
    integer, intent(in) :: nl, ipha
    real(8), intent(in) :: alpha(nl), beta(nl)
    real(8), intent(in) :: rho(nl), h(nl)
    real(8), intent(in) :: omega, p
    logical, intent(in) :: sea_flag
    complex(kind(0d0)), intent(out) :: ur, uz
    complex(kind(0d0)) :: e_inv(4,4), p_mat(4,4), sl(4,4)
    complex(kind(0d0)) :: p_mat2(4,4), lq(2,2)
    complex(kind(0d0)) :: denom, a, b
    integer :: i, i0
    
    

    if (sea_flag) then
       i0 = 2
    else
       i0 = 1
    end if
    
    call e_inverse(omega, rho(nl), alpha(nl), beta(nl), p, e_inv)
    
    call propagator_sol(omega, rho(i0), alpha(i0), beta(i0), p, h(i0), p_mat)
    do i = i0 + 1, nl - 1
       call propagator_sol(omega, rho(i), alpha(i), beta(i), p, h(i), p_mat2)
       p_mat = matmul(p_mat2, p_mat)
    end do
    
    if (sea_flag) then
       call propagator_liq(omega, rho(1), alpha(1), p, h(1), lq)
    end if
    
    sl = matmul(e_inv, p_mat)
    if (.not. sea_flag) then
       denom = sl(3,1) * sl(4,2) - sl(3,2) * sl(4,1)
       if (ipha >= 0) then
          ur = sl(4,2) / denom
          uz = - sl(4,1) / denom
       else 
          ur = - sl(3,2) / denom
          uz = sl(3,1) / denom
       end if
    else
       ! The follwoing equations are derived from boundary condition at 
       ! solid-liquid interface and free surface 
       a = sl(4,2)*lq(1,1) + sl(4,4)*lq(2,1)
       b = sl(3,2)*lq(1,1) + sl(3,4)*lq(2,1)
       if (ipha >= 0) then
          ur = a / (a*sl(3,1) - b*sl(4,1))
          uz = lq(1,1)*sl(4,1) / (b*sl(4,1)-a*sl(3,1))
       else
          ur = - b / (a*sl(3,1) - b*sl(4,1))
          uz = - lq(1,1)*sl(3,1) / (b*sl(4,1)-a*sl(3,1))
       end if
       
    end if
    
    return 
  end subroutine haskell
  !------------------------------------------------------------
  
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
  ! propagator (Aki & Richards, pp. 398, Eq. (3) in Box 9.1)
  !------------------------------------------------------------
  subroutine propagator_sol(omega, rho, alpha, beta, p, z, p_mat)
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
    
    !p_mat(1,1) = 2.d0*beta2*p2*cos_xi + bp*cos_eta
    !p_mat(1,2) = p*( bp/xi*sin_xi - 2.d0*beta2*eta*sin_eta ) * ei
    !p_mat(1,3) = (p2/xi*sin_xi + eta*sin_eta)/(omega*rho)
    !p_mat(1,4) = p*(-cos_xi + cos_eta)/(omega*rho) * ei  
    !p_mat(2,1) = p*( 2.d0*beta2*xi*sin_xi - bp/eta*sin_eta ) * ei
    !p_mat(2,2) = bp*cos_xi + 2.d0*beta2*p2*cos_eta
    !p_mat(2,3) = p_mat(1,4)
    !p_mat(2,4) = (xi*sin_xi + p2/eta*sin_eta)/(omega*rho)
    !p_mat(3,1) = omega*rho*( -4.d0*beta2*beta2*p2*xi*sin_xi - bp*bp/eta*sin_eta )
    !p_mat(3,2) = 2.d0*omega*beta2*rho*p*bp*( cos_xi - cos_eta ) * ei
    !p_mat(3,3) = p_mat(1,1)
    !p_mat(3,4) = p_mat(2,1)
    !p_mat(4,1) = p_mat(3,2)
    !p_mat(4,2) = -omega*rho*( bp*bp/xi*sin_xi + 4.d0*beta2*beta2*p2*eta*sin_eta  )
    !p_mat(4,3) = p_mat(1,2)  
    !p_mat(4,4) = p_mat(2,2)

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
  end subroutine propagator_sol
  
  !------------------------------------------------------------
  subroutine propagator_liq(omega, rho, alpha, p, z, p_mat)
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
  end subroutine propagator_liq
  
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
    
    amp(1:n) =  abs(x(1:n))
    
    wlvl = pcnt * maxval(abs(x(1:n)))
    
    do i = 1, n
       z(i) = y(i) * conjg(x(i))/ max(amp(i), wlvl)
    end do
    
    return
  end subroutine water_level_decon
  
end module forward
