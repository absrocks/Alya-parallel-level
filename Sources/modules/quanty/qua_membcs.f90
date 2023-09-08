subroutine qua_membcs(itask)
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_membcs
  ! NAME
  !    qua_membcs
  ! DESCRIPTION
  !    Allocate memory for the physical problem
  ! OUTPUT 
  ! USES

  ! USED BY
  !    qua_reabcs
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_quanty
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin
  integer(4)              :: istat

  select case(itask)

  case(1_ip)
     !
     ! Fixity and boundary values
     !
     allocate(kfl_fixno_qua(1,npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXNO_QUA','qua_membcs',kfl_fixno_qua)
     do ipoin=1,npoin
        kfl_fixno_qua(1,ipoin)=-1
     end do
     allocate(kfl_fixbo_qua(nboun),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXBO_QUA','qua_membcs',kfl_fixbo_qua)

     if( kfl_conbc_qua == 0 ) then
        allocate(bvess_qua(npoin,2),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_QUA','qua_membcs',bvess_qua)
        allocate(bvessH_qua(npoin,2),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BVESSH_QUA','qua_membcs',bvessH_qua)
     else
        allocate(bvess_qua(npoin,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_QUA','qua_membcs',bvess_qua)
        allocate(bvessH_qua(npoin,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BVESSH_QUA','qua_membcs',bvessH_qua)
     end if

  case(2_ip)
     !
     ! Non-constant b.c.'s : Functions
     !
     allocate(kfl_funno_qua(npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FUNNO_QUA','qua_membcs',kfl_funno_qua)
     allocate(kfl_funbo_qua(nboun),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FUNBO_QUA','qua_membcs',kfl_funbo_qua)

     allocate(kfl_funty_qua(10),   stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FUNTY_QUA','qua_membcs',kfl_funty_qua)
     allocate(funpa_qua(6,10),     stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'FUNPA_QUA',    'qua_membcs',funpa_qua)

  case(3_ip)
     !
     ! Deallocate memory
     !
     call memchk(two,istat,mem_modul(1:2,modul),'KFL_FIXNO_QUA','qua_membcs',kfl_fixno_qua)
     deallocate(kfl_fixno_qua,stat=istat)
     if(istat/=0) call memerr(two,'KFL_FIXNO_TEM','nsi_membcs',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'KFL_FIXBO_QUA','qua_membcs',kfl_fixbo_qua)
     deallocate(kfl_fixbo_qua,stat=istat)
     if(istat/=0) call memerr(two,'KFL_FIXBO_QUA','nsi_membcs',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'BVESS_QUA','qua_membcs',bvess_qua)
     deallocate(bvess_qua,stat=istat)
     if(istat/=0) call memerr(two,'BVESS_QUA','nsi_membcs',0_ip)

  end select

end subroutine qua_membcs
