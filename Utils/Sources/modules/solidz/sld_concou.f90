subroutine sld_concou()
!-----------------------------------------------------------------------
!****f* Solidz/sld_concou
! NAME 
!    sld_concou
! DESCRIPTION
!    This routine checks the SOLIDZ convergence of the run.
! USED BY
!    Solidz
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solidz
  use mod_outfor, only : outfor
  implicit none
  !
  ! Check convergence
  !
  if(kfl_conve(modul)==1) then
     !if(resid_sld>cotol_sld) kfl_gocou = 1
     kfl_gocou = 1
  end if
  glres(modul) = resid_sld
  !
  ! Output residuals
  !
  coutp(1)='DISPLACEMENT'
  routp(1)=resid_sld
  call outfor(9_ip,lun_outpu,' ')
    
end subroutine sld_concou
