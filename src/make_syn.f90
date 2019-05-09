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

program make_syn
  use params
  use model
  use likelihood
  use forward
  use fftw
  use pt_mcmc
  use math
  implicit none 
  integer :: nlay, i, itrc, it, ierr, iarg
  real(8) :: alpha(nlay_max), beta(nlay_max), rho(nlay_max), h(nlay_max)
  logical, parameter :: verb = .true.
  real(8), allocatable :: noise(:,:), noise_sigma(:)
  character(clen_max) :: out_sac, param_file
  logical :: is_valid

  
  param_file = "params.in"
  iarg = command_argument_count()
  if (iarg > 0) then
     call get_command_argument(1, param_file)   
   end if

  call get_params(verb, param_file)
  
  call sgrnd(iseed)

  call read_ref_model(verb)
  
  call init_model(verb)

  call read_obs(verb)

  call init_fftw()

  call init_forward(verb)

  call init_likelihood(verb)
  
    
  call format_model(k(1), z(:,1), dvp(:,1), dvs(:,1), &
       & nlay, alpha, beta, rho, h, is_valid)

  open(54, file = "test_vel", status = "unknown", iostat = ierr)
  if (ierr /= 0) then 
     write(0,*)"ERROR: cannot create test_vel"
     stop
  end if
  do i = 1, nlay
     write(54,*)alpha(i), beta(i), rho(i), h(i)
  end do
  close(54)
  
  ! Add noise
  allocate(noise(nfft, ntrc), noise_sigma(ntrc))

  
  if (is_ray_common) then
     noise_sigma(1) = grnd() * (sig_max(1) - sig_min(1)) + sig_min(1)
     do it = 1, nfft
        noise(it, 1) = gauss() * noise_sigma(1)
     end do
     
     do itrc = 1, ntrc
        rx(1:nfft) = noise(1:nfft, 1)
        call dfftw_execute(ifft2)
        cx(1:nfft/2+1) = cx(1:nfft/2+1) * flt(1:nfft/2+1, itrc)
        call dfftw_execute(ifft)
        noise(1:nfft, itrc) = rx(1:nfft)
     end do
     
     write(*,*)"Noise level of all traces: ", noise_sigma(1)
     
  else
     do itrc = 1, ntrc
        noise_sigma(itrc) = grnd() * &
             &(sig_max(itrc) - sig_min(itrc)) + sig_min(itrc)
        do it = 1, nfft
           noise(it, itrc) = gauss() * noise_sigma(itrc)
        end do
        
        rx(1:nfft) = noise(1:nfft, itrc)
        call dfftw_execute(ifft2)
        cx(1:nfft/2+1) = cx(1:nfft/2+1) * flt(1:nfft/2+1, itrc)
        call dfftw_execute(ifft)
        noise(1:nfft, itrc) = rx(1:nfft)
        
        write(*,*)"Noise level of trace ", itrc, ":", noise_sigma(itrc)
     end do
  end if

  ! output to SAC
  do itrc = 1, ntrc
     write(out_sac,'(A10,I2.2, A2)') "test_trace.",  itrc, "wn"
     open(55, file = out_sac, status = "unknown", access = "direct", &
          & recl = 4, iostat = ierr)
     if (ierr /= 0) then
        write(0,*)"ERROR: cannot create ", trim(out_sac)
        stop
     end if
     write(55, rec = 1) real(delta, kind(0e0))
     write(55, rec = 6) real(t_start, kind(0e0))
     write(55, rec = 7) real(t_end, kind(0e0))
     write(55, rec = 77) 6
     write(55, rec = 86) 1
     write(55, rec = 80) nsmp
     write(55, rec = 106) 1
     do it = 1, nsmp
        write(55, rec = 158 + it) &
             & real(rft(it, itrc, 1) + noise(it, itrc), kind(0e0))
     end do
     close(55)

     write(out_sac,'(A10,I2.2)') "test_trace.",  itrc
     open(56, file = out_sac, status = "unknown", access = "direct", &
          & recl = 4, iostat = ierr)
     if (ierr /= 0) then
        write(0,*)"ERROR: cannot create ", trim(out_sac)
        stop
     end if
     write(56, rec = 1) real(delta, kind(0e0))
     write(56, rec = 6) real(t_start, kind(0e0))
     write(56, rec = 7) real(t_end, kind(0e0))
     write(56, rec = 77) 6
     write(56, rec = 86) 1
     write(56, rec = 80) nsmp
     write(56, rec = 106) 1
     do it = 1, nsmp
        write(56, rec = 158 + it) &
             & real(rft(it, itrc, 1), kind(0e0))
     end do
     close(56)

     
  end do

end program make_syn



