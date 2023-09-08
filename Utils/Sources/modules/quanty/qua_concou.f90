subroutine qua_concou()
!-----------------------------------------------------------------------
!****f* Quanty/qua_concou
! NAME 
!    qua_concou
! DESCRIPTION
!    This routine checks the convergence of the run.
! USED BY
!    Quanty
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_quanty
  use mod_outfor, only : outfor
  implicit none
  !
  ! Check convergence
  !
  if(kfl_conve(modul)==1) then
     if(resid_qua>cotol_qua) kfl_gocou = 1 
  end if
  glres(modul) = resid_qua
  !
  ! Output residuals
  !
  coutp(1)='Eigenstate'
  routp(1)=resid_qua
  call outfor(9_ip,lun_outpu,' ')

end subroutine qua_concou
