subroutine reaous()
  !------------------------------------------------------------------------
  !****f* kernel/reamod
  ! NAME 
  !    reamod
  ! DESCRIPTION
  !    This routine reads postprocess
  ! USES
  ! USED BY
  !    ***_reapro
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: dummi

  if( INOTSLAVE ) then
     !
     ! Reach the section
     !
     call ecoute('reaous')
     do while(words(1)/='OUTPU')
        call ecoute('reaous')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDOU')
        call ecoute('reaous')
        call posdef(2_ip,dummi)
     end do
  end if

end subroutine reaous
