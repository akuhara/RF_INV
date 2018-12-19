module fftw
  use params, only: nsmp
  implicit none 
  include 'fftw3.f'
  
  complex(kind(0d0)), allocatable :: cx(:)
  real(kind(0d0)), allocatable :: rx(:)
  integer(8) :: ifft_r2c, ifft_c2r

  
contains
  
  !=====================================================================
  subroutine init_fftw()
    
    allocate(cx(nsmp), rx(nsmp))
    call dfftw_plan_dft_r2c_1d(ifft_r2c, nsmp, rx, cx, FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_1d(ifft_c2r, nsmp, cx, rx, FFTW_ESTIMATE)
    
    return 
  end subroutine init_fftw
  !=====================================================================

end module fftw
