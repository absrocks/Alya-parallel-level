subroutine qua_memose()
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_memose
  ! NAME 
  !    qua_memose
  ! DESCRIPTION
  !    Estos sets son de postproceso. No me preocupan por ahora
  !    Allocate memory for element, boundary and node sets
  !    The size of element and boundary set (IESET,IBSET) is stored in:
  !    - VESET_QUA(IESET,NVARS_TEM+1)
  !    - VBSET_QUA(IBSET,NVARS_TEM+1)
  ! USES
  ! USED BY
  !    qua_memall
  !***
  !----------------------------------------------------------------------- 
  use def_parame
  use def_master
  use def_domain
  use def_quanty
  use mod_iofile
  use mod_memchk
  implicit none
  integer(4) :: istat

  if(maxval(npp_setse_qua)>0) then
     allocate(veset_qua(nvars_qua+1,neset),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VESET_QUA','qua_memose',veset_qua)
  end if
  if(maxval(npp_setsb_qua)>0) then
     allocate(vbset_qua(nvars_qua+1,nbset),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VBSET_QUA','qua_memose',vbset_qua)
  end if
  if(maxval(npp_setsn_qua)>0) then
     allocate(vnset_qua(nvars_qua,nnset),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VNSET_QUA','qua_memose',vnset_qua)
  end if

end subroutine qua_memose
