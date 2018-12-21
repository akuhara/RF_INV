module forward
  implicit none 
  
  real(8), allocatable, private :: flt(:,:)
  real(8), parameter, private :: pi = 3.1415926535897931
  
contains
  
  !=====================================================================
  
  subroutine init_filter()
    use params, only: a_gus, ntrc, nfft, delta
    implicit none 
    integer :: itrc, i, nh
    real(8) :: df

    nh = nfft / 2 + 1
    df = 1.d0 / (delta * nfft)
    
    allocate(flt(nh, ntrc))
    
    do itrc = 1, ntrc
       do i = 1, nh
          flt(i, itrc) = exp(-(2.d0 * pi * df / (2.d0 * a_gus(itrc)))**2)
       end do
    end do
    


    return 
  end subroutine init_filter
  
  !=====================================================================

  
end module forward
