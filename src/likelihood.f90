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

module likelihood
  implicit none 
  real(8), allocatable, public :: sig(:,:)
  real(8), allocatable, public :: rft(:,:,:)
  real(8), allocatable, public :: log_likelihood(:)

  real(8), allocatable, private :: r_inv(:,:,:)
    
  public init_likelihood, calc_likelihood
  private init_sig, init_r_inv, init_rft
  
contains
  !=====================================================================

  subroutine init_likelihood(verb)
    implicit none
    logical, intent(in) :: verb
    
    call init_sig(verb)
    call init_r_inv(verb)
    call init_rft()
    
    
    return 
  end subroutine init_likelihood
    
  !=====================================================================

  subroutine calc_likelihood(chain_id, itype, prop_k, prop_z, &
       & prop_dvp, prop_dvs, sig, prop_log_likelihood, prop_rft)
    use params
    use forward
    use model
    implicit none
    integer, intent(in) :: prop_k, chain_id
    !logical, intent(in) :: fwd_flag
    integer, intent(in) :: itype
    real(8), intent(in) :: prop_z(k_max-1), prop_dvp(k_max) 
    real(8), intent(in) :: prop_dvs(k_max), sig(ntrc)
    real(8), intent(out) :: prop_log_likelihood
    real(8), intent(out) :: prop_rft(nfft, ntrc)
    integer :: itrc, nlay, it
    real(8) :: alpha(nlay_max), beta(nlay_max), rho(nlay_max), h(nlay_max)
    real(8) :: misfits(nsmp), phi1(nsmp), phi, s
    real(8), parameter :: pi = 3.1415926535897931
    logical :: is_valid
    
    if (itype == itype_sig) then
       prop_rft(1:nfft, 1:ntrc) = rft(1:nfft, 1:ntrc, chain_id) 
    else
       call format_model(prop_k, prop_z, prop_dvp, prop_dvs, &
            & nlay, alpha, beta, rho, h, is_valid)
       
       call calc_rf(chain_id, nlay, nfft, ntrc, rayps, &
            & alpha, beta, rho, h, itype, prop_rft)
    end if
    
    

    prop_log_likelihood = 0.d0
    do itrc = 1, ntrc
       misfits(1:nsmp) = prop_rft(1:nsmp,itrc) - obs(1:nsmp, itrc)
       
       ! Likelihood
       s = sig(itrc)
       phi1 = matmul(misfits, r_inv(:,:,itrc))
       phi = dot_product(phi1, misfits)
       prop_log_likelihood = &
            & prop_log_likelihood - 0.5d0 * phi / (s * s) &
            & - dble(nsmp) * log(s)

    end do
    
    return 
  end subroutine calc_likelihood

  
  !=====================================================================
  
  
  subroutine init_sig(verb)
    use params, only: sig_min, sig_max, nchains, ntrc, sig_mode
    use mt19937
    implicit none 
    logical, intent(in) :: verb
    integer :: ichain, itrc, i
    real(8) :: tmp
    
    allocate(sig(ntrc, nchains))

    if (sig_mode == 1) then
       do ichain = 1, nchains
          do itrc = 1, ntrc
             sig(itrc, ichain) = sig_min + grnd() * (sig_max - sig_min)
          end do
       end do
    else 
       do ichain = 1, nchains
          do itrc = 1, ntrc
             sig(itrc, ichain) = sig_min
          end do
       end do
    end if
       
    
    if (verb) then
       write(*,*)
       write(*,*)"--- Initialize noise sigma ---"
       write(*,*)"  sample output"
       do ichain = 1, nchains
          write(*,*)ichain, sig(1, ichain)
       end do
    end if

    return 
  end subroutine init_sig
  
  !=====================================================================
  
  subroutine init_rft()
    use params, only: nfft, ntrc, nchains, itype_birth
    use model, only: k, z, dvp, dvs
    implicit none 
    integer :: ichain
    

    allocate(rft(nfft, ntrc, nchains))
    allocate(log_likelihood(nchains))
    
    do ichain = 1, nchains
       call calc_likelihood(ichain, itype_birth, k(ichain), z(:,ichain), &
            & dvp(:,ichain), dvs(:,ichain), sig(:,ichain), &
            & log_likelihood(ichain), rft(:,:,ichain))
    end do
    
    

    return 
  end subroutine init_rft
  
  !=====================================================================


  subroutine init_r_inv(verb)
    use params
    implicit none
    logical, intent(in) :: verb  
    integer :: i, j, itrc, ierr, lwork
    real(8), allocatable :: r_mat(:,:), s(:), u(:,:), vt(:,:), &
         & work(:), diag(:,:), test(:,:)
    real(8) :: r, tmpr, dummy(1, 1)
    
    
    allocate(r_inv(nsmp, nsmp, ntrc), r_mat(nsmp, nsmp))
    allocate(s(nsmp), u(nsmp, nsmp), vt(nsmp, nsmp))
    allocate(diag(nsmp, nsmp), test(nsmp, nsmp))
    
    do itrc = 1, ntrc
       r = exp(-a_gus(itrc)**2 * delta**2)
       r_mat(1:nsmp, 1:nsmp) = 0.d0
       do i = 1, nsmp
          r_mat(i, i) = 1.d0
       end do
       do i = 1, nsmp - 1
          tmpr = r ** (i * i)
          if (tmpr < 1.0d-3) exit
          do j = 1, nsmp - 1
             r_mat(j, j+1) = tmpr
          end do
       end do
       do i = 2, nsmp - 1
          do j = 1, i - 1
             r_mat(i, j) = r_mat(j, i)
          end do
       end do
       
       
       if (verb) then ! save r_mat value for confirmation below
          test(1:nsmp, 1:nsmp) = r_mat(1:nsmp, 1:nsmp)
       end if
       
       call dgesvd('A', 'A', nsmp, nsmp, r_mat, nsmp, s, u, nsmp, &
            & vt, nsmp, dummy, -1, ierr)
       lwork = nint(dummy(1,1))
       if (itrc == 1) then
          allocate(work(lwork))
       end if
       
       call dgesvd('A', 'A', nsmp, nsmp, r_mat, nsmp, s, u, nsmp, &
            & vt, nsmp, work, lwork, ierr)
       if (ierr /= 0 .and. verb) then
          write(0,*)"ERROR: in calculating inverse R matrix"
          call mpi_finalize(ierr)
          stop
       end if
       
       diag(1:nsmp, 1:nsmp) = 0.d0
       do i = 1, nsmp
          if (s(i) > 1.0d-5) then
             diag(i, i) = 1.d0 / s(i)
          else
             diag(i, i) = 0.d0
          end if
       end do
       
       r_inv(:, :, itrc) = &
            & matmul(matmul(transpose(vt),diag),transpose(u))
       
       ! confirmation
       if (verb) then
          test = &
               & matmul(r_inv(1:nsmp, 1:nsmp, itrc), test(1:nsmp, 1:nsmp))
          write(*,*)
          write(*,*) "--- Confirmation of R inverse matrix ---"
          write(*,*) " * trace # = ", itrc
          write(*,*) " r         = ", r
          do i = 1, nsmp, 100
             write(*,*) test(1, i)
          end do
       end if
       
       
    end do
    
    return 
  end subroutine init_r_inv
  
  !=====================================================================
  
end module likelihood
