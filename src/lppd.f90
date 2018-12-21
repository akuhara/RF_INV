! This module uses LAPACK library for singular value decomposition
module lppd
  implicit none 
  real(8), allocatable :: r_inv(:,:,:)
  real(8), allocatable :: sig(:,:)
  real(8), allocatable :: log_fac(:)
  
contains
  !=====================================================================

  subroutine init_sig(verb)
    use params, only: sig_min, sig_max, nchains, ntrc, k_max
    use mt19937
    implicit none 
    logical, intent(in) :: verb
    integer :: ichain, itrc, i
    real(8) :: tmp
    
    allocate(sig(ntrc, nchains))
    allocate(log_fac(k_max))

    do ichain = 1, nchains
       do itrc = 1, ntrc
          sig(itrc, ichain) = sig_min + grnd() * (sig_max - sig_min)
       end do
    end do
    
    tmp = 0.d0
    do i = 1, k_max
       tmp = tmp + log(dble(i))
       log_fac(i) = tmp
    end do
    
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

  subroutine calc_log(ichain, lppd)
    use params
    use model
    implicit none
    integer, intent(in) :: ichain
    real(8), intent(out) :: lppd
    integer :: itrc, nlay, ki 
    real(8) :: alpha(nlay_max), beta(nlay_max), rho(nlay_max), h(nlay_max)
    real(8) :: misfits(nsmp), phi1(nsmp), phi, s
    real(8), parameter :: pi = 3.1415926535897931
    
    call format_model(ichain, nlay, alpha, beta, rho, h)
    ki = k(ichain)
    

    lppd = 0.d0
    do itrc = 1, ntrc
       s = sig(itrc, ichain)
    
       ! forward modeling
       misfits(1:nsmp) = 0.d0

       ! Likelihood
       phi1 = matmul(misfits,r_inv(:,:,itrc))
       phi = dot_product(phi1,misfits)
       lppd = lppd - 0.5d0 * phi / (s * s) - dble(nsmp) * log(s)
       
       ! Prior
       

    end do
    

    return 
  end subroutine calc_log

  
  !=====================================================================

    
  subroutine calc_r_inv(verb)
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
       allocate(work(lwork))
       
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
  end subroutine calc_r_inv
  
  !=====================================================================
  
end module lppd
