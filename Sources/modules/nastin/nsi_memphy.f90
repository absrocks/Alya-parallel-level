subroutine nsi_memphy(itask)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_memphy
  ! NAME 
  !    nsi_memphy
  ! DESCRIPTION
  !    This routine allocates memory for physical arrays
  ! USES
  !    ecoute
  !    memchk
  !    runend
  ! USED BY
  !    nsi_reaphy
  !***
  !-----------------------------------------------------------------------
  use def_parame 
  use def_inpout
  use def_master
  use def_nastin
  use def_domain
  use mod_memchk
  use mod_memory
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: icoef,imate
  integer(4)              :: istat

  select case(itask)

  case(1_ip)
     !
     ! List of element materials
     ! 

  case(2_ip)
     !
     ! Properties: Allocate memory
     !

  case(3_ip)
     !
     ! Properties: deallocate memory
     ! DO NOT DEALLOCATE BECAUSE THEY ARE USED BY NSI_OUTINF
     !

  case( 4_ip)
     !
     ! Properties: Allocate memory
     !
     if( INOTMASTER ) then
        call memory_alloca(mem_modul(1:2,modul),'PROPE_NSI','nsi_memphy',prope_nsi,2_ip,npoin)
     end if

  case(-4_ip)
     !
     ! Properties: Deallocate memory
     !
     if( INOTMASTER ) then
        call memory_deallo(mem_modul(1:2,modul),'PROPE_NSI','nsi_memphy',prope_nsi)
     end if

  case( 5_ip)
     !
     ! Material force
     !
     call memory_alloca(mem_modul(1:2,modul),'LFORC_MATERIAL_NSI','nsi_memphy',lforc_material_nsi,nmate)
     call memory_alloca(mem_modul(1:2,modul),'XFORC_MATERIAL_NSI','nsi_memphy',xforc_material_nsi,mforc_material_nsi,nmate)
     call memory_alloca(mem_modul(1:2,modul),'VELTA_NSI','nsi_memphy',velta_nsi,mtabl_nsi, nmate)
     call memory_alloca(mem_modul(1:2,modul),'THRTA_NSI','Nsi_memphy',thrta_nsi,mtabl_nsi, nmate)
     call memory_alloca(mem_modul(1:2,modul),'POWTA_NSI','nsi_memphy',powta_nsi,mtabl_nsi, nmate)
     call memory_alloca(mem_modul(1:2,modul),'VEAVE_NSI','nsi_memphy',veave_nsi,mtabl_nsi, nmate)
     call memory_alloca(mem_modul(1:2,modul),'NTABL_NSI','nsi_memphy',ntabl_nsi, nmate)
     call memory_alloca(mem_modul(1:2,modul),'NTABR_NSI','nsi_memphy',ntabr_nsi, nmate)
     ntabr_nsi = 0
     call memory_alloca(mem_modul(1:2,modul),'RADIU_NSI','nsi_memphy',radiu_nsi,mtabl_nsi, nmate)
     call memory_alloca(mem_modul(1:2,modul),'FORCN_NSI','Nsi_memphy',forcn_nsi,mtabl_nsi, nmate)
     call memory_alloca(mem_modul(1:2,modul),'FORCT_NSI','nsi_memphy',forct_nsi,mtabl_nsi, nmate)
     
  case(-5_ip)
     !
     ! Material force : deallocates structures if master
     !
     if (.false.) then
        call memory_deallo(mem_modul(1:2,modul),'LFORC_MATERIAL_NSI','nsi_memphy',lforc_material_nsi)
        call memory_deallo(mem_modul(1:2,modul),'XFORC_MATERIAL_NSI','nsi_memphy',xforc_material_nsi)
        call memory_deallo(mem_modul(1:2,modul),'VELTA_NSI','nsi_memphy',velta_nsi)
        call memory_deallo(mem_modul(1:2,modul),'THRTA_NSI','nsi_memphy',thrta_nsi)
        call memory_deallo(mem_modul(1:2,modul),'POWTA_NSI','nsi_memphy',powta_nsi)
        call memory_deallo(mem_modul(1:2,modul),'VEAVE_NSI','nsi_memphy',veave_nsi)
        call memory_deallo(mem_modul(1:2,modul),'NTABL_NSI','nsi_memphy',ntabl_nsi)
        call memory_deallo(mem_modul(1:2,modul),'NTABR_NSI','nsi_memphy',ntabr_nsi)
     end if
     
  case( 6_ip)
     !
     ! List for boundary nodes
     !
     call memory_alloca(mem_modul(1:2,modul),'BNTAB_NSI','nsi_memphy',bntab_nsi,nbnod_nsi,3_ip)

  case( 7_ip)
     !
     ! List for boundary nodes
     !
     call memory_alloca(mem_modul(1:2,modul),'BNVAL_NSI','nsi_memphy',bnval_nsi,nbval_nsi*nbtim_nsi,3_ip)
     call memory_alloca(mem_modul(1:2,modul),'IBOUN_NSI','nsi_memphy',iboun_nsi,npoin)

  end select

end subroutine nsi_memphy
