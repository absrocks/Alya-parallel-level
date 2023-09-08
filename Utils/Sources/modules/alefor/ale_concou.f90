subroutine ale_concou
!-----------------------------------------------------------------------
!****f* Alefor/ale_concou
! NAME 
!    ale_concou
! DESCRIPTION
!    This routine checks the ALE formulation convergence of the run.
! USED BY
!    Alefor
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_alefor
  use      def_coupli, only : kfl_gozon
  implicit none

  !!if(kfl_conve(modul)==1) then
  !!   if(resid_ale>cotol_ale) kfl_gocou = 1
  !!end if
  !!glres(modul) = resid_ale
  !!write(lun_outpu,300) resid_ale

! 300 format(  10x,'Residual for the temperature: ',3x,e16.6,' (%)')
    
end subroutine ale_concou
