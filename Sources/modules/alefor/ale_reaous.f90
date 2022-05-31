subroutine ale_reaous
  !-----------------------------------------------------------------------
  !****f* Nastin/ale_reaous
  ! NAME 
  !    ale_reaous
  ! DESCRIPTION
  !    This routine reads the output strategy for the incompressible NS
  !    equations.
  ! USES
  !    ecoute
  ! USED BY
  !    ale_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_alefor
  use def_domain
  use mod_memchk
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: dummi

  if( INOTSLAVE ) then
     !
     ! Reach the section
     !
     call ecoute('ale_reaous')
     do while(words(1)/='OUTPU')
        call ecoute('ale_reaous')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDOU')
        call ecoute('ale_reaous')        
        call posdef(2_ip,dummi) 
     end do

  end if

end subroutine ale_reaous

