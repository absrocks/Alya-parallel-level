subroutine nsa_concou
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_concou
  ! NAME 
  !    nsa_concou
  ! DESCRIPTION
  !    This routine checks the convergence of the run and
  !    set the general convergence flags.
  ! USED BY
  !    Nastal
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_nastal
  use mod_commdom_alya, only: INONE
  use mod_outfor, only : outfor
  implicit none

  !call nsa_coupli(ITASK_CONCOU)

  if(kfl_conve(modul)==1) then
     if(resid_nsa(4)>cotol_nsa) kfl_gocou = 1
  end if
  glres(modul) = resid_nsa(4)

  coutp='Momentum'
  routp(1)=resid_nsa(1)
  call outfor(9_ip,lun_outpu,' ')
  coutp='Continuity'
  routp(1)=resid_nsa(2)
  call outfor(9_ip,lun_outpu,' ')
  if (kfl_foreg_nsa == 0) then
     coutp='Energy'
     routp(1)=resid_nsa(3)
     call outfor(9_ip,lun_outpu,' ')     
  end if
  coutp='Global'
  routp(1)=resid_nsa(4)
  call outfor(9_ip,lun_outpu,' ')

 call nsa_coupli(ITASK_CONCOU)

end subroutine nsa_concou
