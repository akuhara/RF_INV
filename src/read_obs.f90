subroutine read_obs(verb)
  use params
  implicit none 
  logical, intent(in) :: verb
  integer :: itrc, ierr, npts, it1, it2, it
  real :: delta4, t_beg4, tmp4(npts_max)
  character(clen_max) :: ofile


  if (verb) then
     write(*,*)
     write(*,*)"--- Reading observed data ---"
  end if

  
  do itrc = 1, ntrc
     open(io_obs, file = obs_files(itrc), status = "old", &
          & access = "direct", recl = 4, iostat = ierr)
     if (ierr /= 0 ) then
        if (verb) write(*,*) "ERROR: cannot open ", &
             & trim(obs_files(itrc))
        call mpi_finalize(ierr)
        stop
     end if
     
     read(io_obs, rec = 1) delta4
     read(io_obs, rec = 6) t_beg4
     read(io_obs, rec = 80) npts
     it1 = nint((t_start - t_beg4) / delta4) + 1
     it2 = nint((t_end - t_beg4) / delta4) + 1
     nsmp = it2 - it1 + 1
     delta = dble(delta4)


     do it = 1, nsmp
        read(io_obs, rec = 158 + it + it1 - 1) tmp4(it)
     end do
     obs(1:nsmp, itrc) = tmp4(1:nsmp)
     close(io_obs)

     write(ofile, '(A5,I2.2)')"input", itrc
     
     open(io_obs, file = ofile, status = "unknown")
     do it = 1, nsmp
        write(io_obs,*)(it - 1) * delta + t_start, obs(it, itrc)
     end do
     close(io_obs)

     if (verb) then
        write(*,*) "Finish reading ", trim(obs_files(itrc))
        write(*,*) "  output -> ", trim(ofile)
     end if
  end do
  
  return 
end subroutine read_obs
