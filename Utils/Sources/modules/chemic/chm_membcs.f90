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
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case(itask)

  case(1_ip)
     !
     ! Fixity and boundary values
     !
     allocate(kfl_fixno_chm(nclas_chm,npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXNO_CHM','chm_membcs',kfl_fixno_chm)

     allocate(kfl_fixbo_chm(nclas_chm,nboun),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXBO_CHM','chm_membcs',kfl_fixbo_chm)

     allocate(bvess_chm(nclas_chm,npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_CHM','chm_membcs',bvess_chm)

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

  case(3_ip)
     !
     ! Allocate parameters
     !
     allocate(kfl_initi_chm(nspec_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_INITI_CHM','chm_membcs',kfl_initi_chm)

     allocate(kfl_usrbc_chm(nclas_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_USRBC_CHM','chm_membcs',kfl_usrbc_chm)

     allocate(xinit_chm(nspec_chm,2),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'XINIT_CHM','chm_membcs',xinit_chm)

     allocate(panat_chm(2,nclas_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'PANAT_CHM','chm_membcs',panat_chm)

  case(-3_ip)
     !
     ! Deallocate parameters
     !
     call memchk(two,istat,mem_modul(1:2,modul),'PANAT_CHM','chm_membcs',panat_chm)
     deallocate(panat_chm,stat=istat)
     if(istat/=0) call memerr(two,'PANAT_CHM','nsi_membcs',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'XINIT_CHM','chm_membcs',xinit_chm)
     deallocate(xinit_chm,stat=istat)
     if(istat/=0) call memerr(two,'XINIT_CHM','nsi_membcs',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'KFL_USRBC_CHM','chm_membcs',kfl_usrbc_chm)
     deallocate(kfl_usrbc_chm,stat=istat)
     if(istat/=0) call memerr(two,'KFL_USRBC_CHM','nsi_membcs',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'KFL_INITI_CHM','chm_membcs',kfl_initi_chm)
     deallocate(kfl_initi_chm,stat=istat)
     if(istat/=0) call memerr(two,'KFL_INITI_CHM','nsi_membcs',0_ip)

  end select

end subroutine chm_membcs
