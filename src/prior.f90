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

module prior
  implicit none 
  
  public log_prior_ratio, laplace
  
contains

  !---------------------------------------------------------------------
  real(8) function log_prior_ratio(x_new, x_old, dev, prior_mode)
    implicit none
    real(8), intent(in) :: x_new, x_old, dev
    integer, intent(in) :: prior_mode
    
    ! Calculate log[p(m')] - log[p(m)] for fixed-D proposals

    if (prior_mode == 1) then
       ! Laplace distribution
       log_prior_ratio = -(abs(x_new) - abs(x_old)) / dev
    else if (prior_mode == 2) then
       ! Gaussian distribution
       log_prior_ratio = -((x_new * x_new) - (x_old * x_old)) &
            & / (2.d0 * dev * dev)
    end if
    
    return 
  end function log_prior_ratio
  !---------------------------------------------------------------------
  

  real(8) function laplace()
    ! Draw random values from laplace distriubtion
    use mt19937
    real(8) :: u1, u1p, u1pp, u1ppp, u2
    real(8) :: a, w
    real(8), parameter :: d = 0.69314718055994529d0 ! ln(2)
    integer :: i_sign, k
    
    ! Step 1
    u1 = grnd()
    u1p = 2.d0 * u1
    
    ! Step 2
    if (u1p < 1.d0) then
       i_sign = 1
       u1pp = 1.d0 - u1p
    else
       i_sign = -1
       u1pp = 2.d0 - u1p
    end if

    ! Step 3
    a = 0.d0
    
    step4_5: do 
       ! Step 4
       u1ppp = 2.d0 * u1pp
       
       ! Step 5
       if (u1ppp >= 1.d0) then
          u1 = u1ppp - 1.d0
          exit step4_5
       else
          a = a + d
          u1pp = u1ppp
       end if
    end do step4_5
    
    step6_9: do
       ! Step 6
       w = d * u1
       laplace = i_sign * (a + w)
       k = 1
       
       step7_8: do
          ! Step 7
          u2 = grnd()
          
          ! Step 8
          if (u2 >= w) then
             u1 = (u2 - w) / (1.d0 - w)
             exit step7_8
          else
             w = u2
             k = k + 1
          end if
       end do step7_8
       
       ! Step 9
       if (mod(k,2) == 1) then
          exit step6_9
       end if
    end do step6_9
    
    ! Step 10
    return 
  end function laplace
  
  

end module prior
