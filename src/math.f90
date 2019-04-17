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

module math
  implicit none 
  
contains
  
  !---------------------------------------------------------------------
  real(8) function gauss()
    use mt19937
    implicit none 
    real(8), parameter :: pi2 = 2.d0 * 3.1415926535897931d0
    real(8) :: v1, v2
    
    v1 = grnd() 
    v2 = grnd()
    gauss = sqrt(-2.d0 * log(v1)) * cos(pi2 * v2)
    
    return 
  end function gauss
  
  !---------------------------------------------------------------------

  
end module math
