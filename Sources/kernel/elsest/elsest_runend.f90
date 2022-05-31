subroutine elsest_runend(message)
  !-----------------------------------------------------------------------
  !
  ! This routine stops the run and writes the summary of CPU time.
  !
  !-----------------------------------------------------------------------
  use def_elsest
  implicit none
  character(*) :: message
  !
  ! Write message and stop the run
  !
  if(iunit(1)/=0) then
!*OMP CRITICAL(write)
     write(iunit(1),1) 'ELSEST ERROR WAS FOUND: '//trim(message)
!*OMP END CRITICAL(write)
  else
!*OMP CRITICAL(write)
     write(6,1)     'ELSEST ERROR WAS FOUND: '//trim(message)
!*OMP END CRITICAL(write)
  end if

  stop

1 format(//,5x,'>>> ',a)

end subroutine elsest_runend
