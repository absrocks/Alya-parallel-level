subroutine rad_concou()
!-----------------------------------------------------------------------
!****f* Radiat/rad_concou
! NAME 
!    rad_concou
! DESCRIPTION
!    This routine checks the radiation convergence of the run.
! USED BY
!    Radiat
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_radiat
  use mod_outfor, only : outfor
  implicit none
  !
  ! Check convergence
  !
  if(kfl_conve(modul)==1) then
     if(resid_rad>cotol_rad) kfl_gocou = 1
  end if
  glres(modul) = resid_rad
  !
  ! Output residuals
  !
  coutp(1)='RADIATION'
  routp(1)=resid_rad
  call outfor(9_ip,lun_outpu,' ')

end subroutine rad_concou
