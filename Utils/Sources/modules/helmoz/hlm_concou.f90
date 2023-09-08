subroutine hlm_concou()
!-----------------------------------------------------------------------
!****f* Helmoz/hlm_concou
! NAME 
!    hlm_concou
! DESCRIPTION
!    This routine checks the helmozature convergence of the run.
! USED BY
!    Helmoz
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_helmoz
  use mod_outfor, only : outfor
  implicit none
  !
  ! Check convergence
  !
  if(kfl_conve(modul)==1) then
     if(resid_hlm>cotol_hlm) kfl_gocou = 1
  end if
  glres(modul) = resid_hlm
  !
  ! Output residuals
  !
  coutp(1)='POTENTIAL'
  routp(1)=resid_hlm
  call outfor(9_ip,lun_outpu,' ')

end subroutine hlm_concou
