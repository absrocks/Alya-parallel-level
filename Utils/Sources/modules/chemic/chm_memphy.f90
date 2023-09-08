subroutine chm_memphy(itask)
  !-----------------------------------------------------------------------
  !****f* partis/chm_memphy
  ! NAME
  !    chm_memphy
  ! DESCRIPTION
  !    Allocate memory for the physical problem
  ! OUTPUT 
  ! USES
  ! USED BY
  !    chm_reaphy
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_chemic
  use mod_memchk
  use def_domain, only : npoin, nelem, ltype, ngaus
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat
  integer(ip)             :: ireac, ielem, pelty, pgaus

  select case(itask)

  case(1_ip)
     !
     ! Allocate properties
     !
     allocate(diffu_chm(2,nclas_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DIFFU_PTS','chm_memphy',diffu_chm)
     allocate(lawdi_chm(2,nclas_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LAWDI_PTS','chm_memphy',lawdi_chm)
     allocate(interaction_chm(nclas_chm,nclas_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'INTERACTION_PTS','chm_memphy',interaction_chm)

     !
     ! Only useful for DEFECT EVOLUTION model
     !
     allocate(equil_chm(2,nclas_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'EQUIL_CHM','chm_memphy',equil_chm)
     !
     ! Only useful for METEO model
     !
     allocate(diame_chm(nclas_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DIAME_CHM','pts_memphy',diame_chm)
     allocate(rhopa_chm(nclas_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'RHOPA_CHM','pts_memphy',rhopa_chm)
     allocate(shape_chm(nclas_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'SHAPE_CHM','pts_memphy',shape_chm)
     allocate(spher_chm(nclas_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'SPHER_CHM','pts_memphy',spher_chm)
     allocate(fract_chm(nclas_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'FRACT_CHM','pts_memphy',fract_chm)

     allocate(speci(nclas_chm)) !!**
     !!FER: CANNOT CHECK ALLOCATION OF THE NEWLY DEFINED TYPE, IT IS NOT IN MEMCHK
     !!!FER WHY  call memchk(zero,istat,mem_modul(1:2,modul),'SPECI_CHM','chm_memphy',speci_chm)

  case(2_ip)
     !
     ! Deallocate properties
     !
     call memchk(two,istat,mem_modul(1:2,modul),'diffu_chm','chm_memphy',diffu_chm)
     deallocate(diffu_chm,stat=istat)
     if(istat/=0) call memerr(two,'diffu_chm','chm_memphy',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'lawdi_chm','chm_memphy',lawdi_chm)
     deallocate(lawdi_chm,stat=istat)
     if(istat/=0) call memerr(two,'lawdi_chm','chm_memphy',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'equil_chm','chm_memphy',equil_chm)
     deallocate(equil_chm,stat=istat)
     if(istat/=0) call memerr(two,'equil_chm','chm_memphy',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'radiu_chm','chm_memphy',radiu_chm)
     deallocate(radiu_chm,stat=istat)
     if(istat/=0) call memerr(two,'radiu_chm','chm_memphy',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'fract_chm','chm_memphy',fract_chm)
     deallocate(fract_chm,stat=istat)
     if(istat/=0) call memerr(two,'fract_chm','chm_memphy',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'spher_chm','chm_memphy',spher_chm)
     deallocate(spher_chm,stat=istat)
     if(istat/=0) call memerr(two,'spher_chm','chm_memphy',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'shape_chm','chm_memphy',shape_chm)
     deallocate(shape_chm,stat=istat)
     if(istat/=0) call memerr(two,'shape_chm','chm_memphy',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'rhopa_chm','chm_memphy',rhopa_chm)
     deallocate(rhopa_chm,stat=istat)
     if(istat/=0) call memerr(two,'rhopa_chm','chm_memphy',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'diame_chm','chm_memphy',diame_chm)
     deallocate(diame_chm,stat=istat)
     if(istat/=0) call memerr(two,'diame_chm','chm_memphy',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'interaction_chm','chm_memphy',interaction_chm)
     deallocate(interaction_chm,stat=istat)
     if(istat/=0) call memerr(two,'interaction_chm','chm_memphy',0_ip)


     call memchk(two,istat,mem_modul(1:2,modul),'effic_chm','chm_memphy',effic_chm)
     deallocate(effic_chm,stat=istat)
     if(istat/=0) call memerr(two,'effic_chm','chm_memphy',0_ip)
    !!!FER WHY    call memchk(two,istat,mem_modul(1:2,modul),'speci_chm','chm_memphy',speci_chm)
     deallocate(speci,stat=istat) !!**
     if(istat/=0) call memerr(two,'speci','chm_memphy',0_ip)

  case(3_ip)
     !
     ! Allocate reaction
     !
     allocate(lreac_chm(nreac_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_memphy',lreac_chm)
     do ireac  = 1,nreac_chm
        allocate(lreac_chm(ireac)%l(ncoef_chm),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_memphy',lreac_chm(ireac)%l)       
     enddo
     allocate(react_chm(ncoef_chm,nreac_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'REACT_CHM','chm_memphy',react_chm)

     allocate(order_chm(nclas_chm,nreac_chm,2),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ORDER_CHM','chm_memphy',order_chm)

     allocate(stoic_chm(nclas_chm,nreac_chm,2),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'STOIC_CHM','chm_memphy',stoic_chm)

     allocate(effic_chm(nclas_chm,nreac_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'EFFIC_CHM','chm_memphy',effic_chm)
  
     allocate(radiu_chm(nreac_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'RADIU_CHM','chm_memphy',radiu_chm)

  case(4_ip)
     !
     ! Deallocate reaction
     !
     call memchk(two,istat,mem_modul(1:2,modul),'lreac_chm','chm_memphy',lreac_chm)
     deallocate(lreac_chm,stat=istat)
     if(istat/=0) call memerr(two,'lreac_chm','chm_memphy',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'react_chm','chm_memphy',react_chm)
     deallocate(react_chm,stat=istat)
     if(istat/=0) call memerr(two,'react_chm','chm_memphy',0_ip)

  case(5_ip)
     !
     ! CFI model: Allocate table
     !
     allocate(table_cfi(1)%ivcfi(table_cfi(1)%ndcfi,maxval(table_cfi(1)%nvcfi)),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'table_cfi','chm_memphy',table_cfi(1)%ivcfi)

     allocate(table_cfi(1)%table(table_cfi(1)%nrcfi,table_cfi(1)%nccfi),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'table_cfi','chm_memphy',table_cfi(1)%table)

     allocate(table_cfi(1)%ymass(table_cfi(1)%nvcfi(3),3),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'table_cfi','chm_memphy',table_cfi(1)%ymass)

  case(6_ip)
     !
     ! CFI model: Allocate table structure
     !
     allocate(table_cfi(1),stat=istat)

  end select

end subroutine chm_memphy
 
