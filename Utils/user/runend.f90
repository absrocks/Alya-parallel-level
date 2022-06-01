subroutine runend(message)
  !-----------------------------------------------------------------------
  !
  ! This routine stops the run and writes the summary of CPU time.
  !
  !-----------------------------------------------------------------------

  implicit none
  character(*) :: message

  
  if( message /= 'O.K.!' ) then 
     print*,'RUNEND MESSAGE|'
     print*,'------>',message
     flush(6)
  end if


end subroutine runend
