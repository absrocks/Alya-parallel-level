subroutine runend(message) 
  !-----------------------------------------------------------------------
  !
  ! This routine stops the run and writes the summary of CPU time.
  !
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use mod_messages
  use mod_outfor,   only : outfor
  use mod_messages, only : livinf
  use mod_messages, only : messages_live
  use mod_iofile,   only : iofile_flush_unit

  implicit none
  character(*)              :: message
  character(1000)           :: my_message
  logical(lg)               :: lopen
  integer                   :: cy

#ifdef makefail 
  integer(ip)               :: kfail,kk
  character(3)              :: tomakefail
     !
     ! To get the line where it has failed 
     !
     kfail=index(message,'O.K')
     if(kfail == 0) then
        tomakefail=' -1'
        read(tomakefail,*) kfail
        write(kfl_paral+6100,*)'message',message
        call iofile_flush_unit(kfl_paral+6100)
        kk = 10_ip / (1_ip+kfail)
        kfail = log (dble(kfail))  ! this one does not fail if there are no debug option the upper one does
     end if
#endif
  !
  ! Write message and stop the run
  !
  if( INOTSLAVE .or. kfl_outpu_par == 1 ) then

     if( message /= 'O.K.!' ) then
        inquire(unit=lun_outpu,opened=lopen)
        if(lopen) call outfor(0_ip,lun_outpu,trim(message))
        if( INOTSLAVE ) then
           call messages_live('|---->  ABORTED IN '//trim(message))
           !
           ! Live a little, a haiku about Alya
           !
           ! Fer
           ! 
           call system_clock(cy) 
           if (mod(cy,100) == 0 ) then
              call livinf(-1_ip,'.',zero)
              call livinf(-1_ip,'..',zero)
              call livinf(-1_ip,'   Yesterday it worked.',zero)
              call livinf(-1_ip,'   Tomorrow too, today not.',zero)
              call livinf(-1_ip,'   Alya is like that.',zero)
              call livinf(-1_ip,'..',zero)
              call livinf(-1_ip,'.',zero)
           endif
        endif
     else
        call outfor(14_ip,lun_outpu,' ')
        call messages_live('')
        call messages_live('CALCULATIONS CORRECT')
     end if

  else if( ISLAVE .or. kfl_outpu_par == 1 ) then

     if( message /= 'O.K.!' ) then
        print*,'RUNEND MESSAGE|'
        print*,'------>',message
        call iofile_flush_unit(6_ip)
     end if

  end if
  !
  ! Personnalized message
  !
  if( INOTSLAVE ) then
     call messages_general(my_message)
     if( trim(my_message) /= '' ) then
        call livinf(-1_ip,'...',zero)
        call livinf(-1_ip,trim(my_message),zero)
        call livinf(-1_ip,'...',zero)
     end if
  end if
  
  if( IPARALL ) then
     if( message == 'O.K.!' ) then
        call par_finali(1_ip) 
     else
        call par_finali(0_ip) 
     end if
  end if

  if( kfl_stop /= 0 ) then
     stop 1
  else
     stop
  end if

end subroutine runend
