subroutine chm_outerr()
!------------------------------------------------------------------------
!****f* partis/chm_outerr
! NAME 
!    chm_outerr
! DESCRIPTION
!    This routine checks if there are errros and warnings
! USES
! USED BY
!    chm_turnon
!***
!------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_chemic
  use mod_outfor, only : outfor
  implicit none
  integer(ip)   :: ierro=0,iwarn=0
  integer(ip)   :: ispec,ireac
  logical       :: l_problem
  character(20) :: messa
  !
  ! TRANSIENT PROBLEM 
  !
  if( kfl_timei /= 0 .and. kfl_timei_chm == 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,'STEADY PARTIS IN A TRANSIENT CALCULATION')
  end if
  !
  ! TEMPERATURE FUNCTION
  !
  if( lawte_chm /= 0 .and. kfl_assem_chm == 3 ) then
     kfl_assem_chm = 2
     iwarn = iwarn+1
     call outfor(2_ip,momod(modul)%lun_outpu,'CANNOT SAVE MATRICES WHEN USING TEMPERATURE FUNCTION')
  end if
  !
  ! LHS Assembly of Source terms: Unable to treat negative Arrhenius powers 
  !
  if (kfl_model_chm == 4_ip .and. kfl_lhsas_chm /= 0 ) then
     l_problem = .false.
     do ireac=1,nreac_chm
        do ispec=1,nspec_chm
           l_problem = l_problem .or. (order_chm(ispec,ireac,1) /= 0.0_rp .and. order_chm(ispec,ireac,1) < 1.0)
           l_problem = l_problem .or. (order_chm(ispec,ireac,2) /= 0.0_rp .and. order_chm(ispec,ireac,2) < 1.0)
        enddo
     enddo
     if( l_problem ) then
        ierro = ierro + 1
        call outfor(2_ip,momod(modul)%lun_outpu,'LHS Assembly of source terms requires Arrhenius exponents >= 1')
     end if
  endif

  !----------------------------------------------------------------------
  !
  ! ERROR MESSAGE
  !
  !----------------------------------------------------------------------
  
  call errors(3_ip,ierro,iwarn,'NULL')

end subroutine chm_outerr
