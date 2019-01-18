subroutine fwd_seis (nlay, npts, rayp, ipha, &
     & alpha, beta, rho, h, ur_freq, uz_freq)
  use params, only: pi2, delta, ei
  implicit none
  integer, intent(in) :: nlay, ipha, npts
  real(8), intent(in) :: alpha(nlay), beta(nlay)
  real(8), intent(in) :: rho(nlay), h(nlay)
  real(8), intent(in) :: rayp
  complex(kind(0d0)), intent(out) :: ur_freq(npts), uz_freq(npts)
  real(8) :: omg
  integer :: i, nhalf
  logical :: sea_flag

  if (beta(1) < 0) then
     sea_flag = .true.
  else
     sea_flag = .false.
  end if
  nhalf = npts / 2 + 1
  
  ! main part by propagator matrix
  !ur_freq(1) = (0.d0, 0.d0)
  !uz_freq(1) = (0.d0, 0.d0)
  do i = 1, nhalf
     omg = pi2 * dble(i-1) / (npts * delta)
     if (omg < 1.0e-5) then
        omg = 1.0e-5
     end if
     call haskell(omg, ur_freq(i), uz_freq(i), rayp, nlay, ipha, &
          & alpha, beta, rho, h, sea_flag)
  end do
  
  return 
end subroutine fwd_seis
    


!-----------------------------------------------------------------------
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
  use params, only: ei
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
  use params, only: ei
  implicit none
  real(8) :: alpha, beta
  real(8), intent(in) :: omega, rho, p, z
  complex(kind(0d0)), intent(out) :: p_mat(4,4)
  real(8) :: eta,xi,beta2,p2,bp,cos_xi,cos_eta,sin_xi,sin_eta
  
  p_mat(:,:) = (0.d0, 0.d0)
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
  p_mat(1,2) = p*( bp/xi*sin_xi - 2.d0*beta2*eta*sin_eta ) * ei
  p_mat(1,3) = (p2/xi*sin_xi + eta*sin_eta)/(omega*rho)
  p_mat(1,4) = p*(-cos_xi + cos_eta)/(omega*rho) * ei  
  p_mat(2,1) = p*( 2.d0*beta2*xi*sin_xi - bp/eta*sin_eta ) * ei
  p_mat(2,2) = bp*cos_xi + 2.d0*beta2*p2*cos_eta
  p_mat(2,3) = p_mat(1,4)
  p_mat(2,4) = (xi*sin_xi + p2/eta*sin_eta)/(omega*rho)
  p_mat(3,1) = omega*rho*( -4.d0*beta2*beta2*p2*xi*sin_xi - bp*bp/eta*sin_eta )
  p_mat(3,2) = 2.d0*omega*beta2*rho*p*bp*( cos_xi - cos_eta ) * ei
  p_mat(3,3) = p_mat(1,1)
  p_mat(3,4) = p_mat(2,1)
  p_mat(4,1) = p_mat(3,2)
  p_mat(4,2) = -omega*rho*( bp*bp/xi*sin_xi + 4.d0*beta2*beta2*p2*eta*sin_eta  )
  p_mat(4,3) = p_mat(1,2)  
  p_mat(4,4) = p_mat(2,2)

  return 
end subroutine propagator_sol

!------------------------------------------------------------
subroutine propagator_liq(omega, rho, alpha, p, z, p_mat)
  use params, only: ei
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


