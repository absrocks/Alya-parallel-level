subroutine chm_membcs(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/chm_membcs
  ! NAME
  !    chm_membcs
  ! DESCRIPTION
  !    Allocate memory for the physical problem
  ! OUTPUT 
  ! USES
  ! USED BY
  !    chm_reaphy
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_chemic
  use mod_memchk
  use mod_memory,         only : memory_alloca
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case(itask)

  case(1_ip)
     !
     ! Fixity and boundary values
     !
     if (kfl_model_chm == 4) then
        allocate(kfl_fixno_chm(nvar_CMC_chm,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXNO_CHM','chm_membcs',kfl_fixno_chm)

        allocate(kfl_fixbo_chm(nvar_CMC_chm,nboun),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXBO_CHM','chm_membcs',kfl_fixbo_chm)

        if (kfl_rstar == 0 .and. kfl_weigh_in_eq_CMC_chm == 1) then
           allocate(bvess_chm(nvar_CMC_chm+1,npoin),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_CHM','chm_membcs',bvess_chm)
        else
           allocate(bvess_chm(nvar_CMC_chm,npoin),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_CHM','chm_membcs',bvess_chm)
        end if

        allocate(bvess_CMC_chm(nZ_CMC_chm,npoin,nvar_CMC_chm),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_CMC_CHM','chm_membcs',bvess_CMC_chm)

        allocate(bvess_ufield_CMC_chm(npoin,nvar_CMC_chm),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_UNCOND_FIELD_CMC_CHM','chm_membcs',bvess_ufield_CMC_chm)

     else
        allocate(kfl_fixno_chm(nclas_chm,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXNO_CHM','chm_membcs',kfl_fixno_chm)
   
        allocate(kfl_fixbo_chm(nclas_chm,nboun),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXBO_CHM','chm_membcs',kfl_fixbo_chm)

        allocate(bvess_chm(nclas_chm,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_CHM','chm_membcs',bvess_chm)
     end if

     kfl_fixno_chm(:,:)           = 0.0_rp
     kfl_fixbo_chm(:,:)           = 0.0_rp
     bvess_chm(:,:)               = 0.0_rp
     if (kfl_model_chm == 4) then
        bvess_CMC_chm(:,:,:)      = 0.0_rp
        bvess_ufield_CMC_chm(:,:) = 0.0_rp
     end if

  case(2_ip)
     !
     ! Deallocate memory
     !
     call memchk(two,istat,mem_modul(1:2,modul),'KFL_FIXNO_CHM','chm_membcs',kfl_fixno_chm)
     deallocate(kfl_fixno_chm,stat=istat)
     if(istat/=0) call memerr(two,'KFL_FIXNO_CHM','nsi_membcs',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'KFL_FIXBO_CHM','chm_membcs',kfl_fixbo_chm)
     deallocate(kfl_fixbo_chm,stat=istat)
     if(istat/=0) call memerr(two,'KFL_FIXBO_CHM','nsi_membcs',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'BVESS_CHM','chm_membcs',bvess_chm)
     deallocate(bvess_chm,stat=istat)
     if(istat/=0) call memerr(two,'BVESS_CHM','nsi_membcs',0_ip)

     if (kfl_model_chm == 4) then
        call memchk(two,istat,mem_modul(1:2,modul),'BVESS_CMC_CHM','chm_membcs',bvess_CMC_chm)
        deallocate(bvess_CMC_chm,stat=istat)
        if(istat/=0) call memerr(two,'BVESS_CMC_CHM','nsi_membcs',0_ip)

        call memchk(two,istat,mem_modul(1:2,modul),'BVESS_UNCOND_FIELD_CMC_CHM','chm_membcs',bvess_ufield_CMC_chm)
        deallocate(bvess_ufield_CMC_chm,stat=istat)
        if(istat/=0) call memerr(two,'BVESS_UNCOND_FIELD_CMC_CHM','nsi_membcs',0_ip)
     end if

  case(3_ip)
     !
     ! Allocate parameters
     !
     allocate(kfl_initi_chm(nclas_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_INITI_CHM','chm_membcs',kfl_initi_chm)

     allocate(xinit_chm(nclas_chm,2),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'XINIT_CHM','chm_membcs',xinit_chm)

  case(-3_ip)
     !
     ! Deallocate parameters
     !
     call memchk(two,istat,mem_modul(1:2,modul),'XINIT_CHM','chm_membcs',xinit_chm)
     deallocate(xinit_chm,stat=istat)
     if(istat/=0) call memerr(two,'XINIT_CHM','nsi_membcs',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'KFL_INITI_CHM','chm_membcs',kfl_initi_chm)
     deallocate(kfl_initi_chm,stat=istat)
     if(istat/=0) call memerr(two,'KFL_INITI_CHM','nsi_membcs',0_ip)
  
  case(4_ip)

     !
     ! Non-constant b.c.'s : Functions
     !
     if (kfl_model_chm == 4) then
        call memory_alloca(mem_modul(1:2,modul),'KFL_FUNNO_CHM','chm_membcs',kfl_funno_chm,npoin,nvar_CMC_chm)
        call memory_alloca(mem_modul(1:2,modul),'KFL_FUNTN_CHM','chm_membcs',kfl_funtn_chm,npoin,nvar_CMC_chm)
     else
        call memory_alloca(mem_modul(1:2,modul),'KFL_FUNNO_CHM','chm_membcs',kfl_funno_chm,npoin,nclas_chm)
        call memory_alloca(mem_modul(1:2,modul),'KFL_FUNTN_CHM','chm_membcs',kfl_funtn_chm,npoin,nclas_chm)
     end if

  end select

end subroutine chm_membcs
