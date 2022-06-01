subroutine runend(message)
  !-----------------------------------------------------------------------
  !
  ! This routine stops the run and writes the summary of CPU time.
  !
  !-----------------------------------------------------------------------
  use      def_kintyp
  implicit none
  character(*) :: message
  logical(lg)  :: lopen
  !
  ! Write message and stop the run
  !
  print*,'RUNEND MESSAGE:'
  print*,message

end subroutine runend
