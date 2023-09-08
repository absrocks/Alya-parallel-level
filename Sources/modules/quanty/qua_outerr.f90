subroutine qua_outerr()
!------------------------------------------------------------------------
!****f* Quanty/qua_outerr
! NAME 
!    qua_outerr
! DESCRIPTION
!    This routine checks if there are errros and warnings
! USES
! USED BY
!    qua_turnon
!***
!------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_quanty
  use mod_outfor, only : outfor
  implicit none
  integer(ip)   :: ierro=0,iwarn=0
  character(20) :: messa


  messa=intost(ierro)
  if(ierro/=0) call outfor(4_ip,lun_outpu_qua,trim(messa))

end subroutine qua_outerr
