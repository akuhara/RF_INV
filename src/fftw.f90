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

module fftw
  use params, only: nfft
  implicit none 
  include 'fftw3.f'
  
  complex(kind(0d0)), allocatable :: cx(:)
  real(kind(0d0)), allocatable :: rx(:)
  integer(8) :: ifft, ifft2

  
contains
  
  !=====================================================================
  subroutine init_fftw()

    allocate(cx(nfft), rx(nfft))
    call dfftw_plan_dft_c2r_1d(ifft, nfft, cx, rx, FFTW_ESTIMATE)
    call dfftw_plan_dft_r2c_1d(ifft2, nfft, rx, cx, FFTW_ESTIMATE)
    
    return 
  end subroutine init_fftw
  !=====================================================================

end module fftw
