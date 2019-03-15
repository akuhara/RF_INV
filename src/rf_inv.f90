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
   use mt19937
   use fftw
   use model
   use likelihood
   use forward
   use pt_mcmc
   use mcmc_out
   implicit none 
   include "mpif.h"
   integer :: nproc, rank, ierr, iarg
   logical :: verb
   character(clen_max) :: param_file

   ! Initialize MPI 
   call mpi_init(ierr)
   call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
   call mpi_comm_rank(MPI_COMM_WORLD, rank,  ierr)
   !write(*,*)nproc, rank
   !if (nproc < 2) then
   !   write(0,*)"ERROR: at least 2 processor is necesarry!"
   !   call mpi_finalize(ierr)
   !   stop
   !end if

   ! Set verbose mode for rank 0
   verb = .false.
   if (rank == 0) verb = .true.  

   !============================================================
   ! Initialize
   !============================================================  
   ! Read parameters from file

   param_file = "params.in"
   iarg = command_argument_count()
   if (iarg > 0) then
      call get_command_argument(1, param_file)   
   end if
   call get_params(verb, param_file) 
   
   ! Read observed files
   call read_obs(verb)

   ! Initialize random number generator
   iseed = iseed + rank * rank * 10000 + 23 * rank
   call sgrnd(iseed)

   ! Initialize FFTW
   call init_fftw()

   call init_forward(verb)
   ! Read reference velocity model
   call read_ref_model(verb)

   ! Generate initial model
   call init_model(verb)
   
   call init_likelihood(verb)
   
   ! Initialize temperature
   call init_pt_mcmc(verb)
   
   !============================================================
   ! MCMC
   !============================================================  

   ! MCMC samping
   call pt_control(verb)
   

   !output
   call output_results(nproc, rank, verb)
   
   ! Finish
   call mpi_finalize(ierr)

   stop
 end program main
