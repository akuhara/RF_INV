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
program main
  use params
  use forward
  use fftw
  implicit none 

  integer, parameter :: nlay = 2
  real(8) :: alpha(nlay), beta(nlay), rho(nlay), h(nlay)
  real(8) :: ray_parameters(1)
  real(8), allocatable :: rft(:)
  integer :: i
  
  sdep = 1.0d0 !
  bdep = 0.0d0
  ntrc = 1
  allocate(a_gus(ntrc), ipha(ntrc))
  a_gus = [8.0d0]
  ipha = [1]
  deconv_mode = 0
  delta = 0.05d0
  t_start = -10.d0
  nfft = 1024
  alpha = [1.5, 5.0]
  beta  = [-9.0,  2.5]
  rho   = [1.0, 3.0]
  h     = [3.0, -10.0]
  ray_parameters = [0.06]
  allocate(rft(nfft))
  
  call init_fftw()

  call init_forward(verb=.true.)

  call calc_rf(1, nlay, nfft, 1, ray_parameters, alpha, beta, rho, h, rft)

  do i = 1, nfft
     write(111,*)(i-1) * delta + t_start, rft(i)
  end do
  stop
end program main
