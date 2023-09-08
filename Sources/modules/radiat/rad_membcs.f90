subroutine rad_membcs(itask)
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_membcs
  ! NAME
  !    rad_membcs
  ! DESCRIPTION
  !    Allocate memory for the physical problem
  ! OUTPUT 
  ! USES
  ! USED BY
  !    rad_reaphy
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_radiat
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin
  integer(4)              :: istat

  select case(itask)

  case(1)
     !
     ! Fixity and boundary values
     !
     allocate(kfl_fixno_rad(1,npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXNO_RAD','rad_membcs',kfl_fixno_rad)
     do ipoin=1,npoin
        kfl_fixno_rad(1,ipoin)=-1
     end do
     allocate(kfl_fixbo_rad(nboun),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXBO_RAD','rad_membcs',kfl_fixbo_rad)
     if(kfl_conbc_rad==0) then
        allocate(bvess_rad(npoin,2),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_RAD','rad_membcs',bvess_rad)
        if(nboun/=0) then
           allocate(bvnat_rad(npnat_rad,nboun,2),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'BVNAT_RAD','rad_membcs',bvnat_rad)
        end if
     else
        allocate(bvess_rad(npoin,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_RAD','rad_membcs',bvess_rad)
        if(nboun/=0) then
           allocate(bvnat_rad(npnat_rad,nboun,1),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'BVNAT_RAD','rad_membcs',bvnat_rad)
        end if
     end if

  case(2)
     !
     ! Non-constant b.c.'s : Functions
     !
     allocate(kfl_funno_rad(npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FUNNO_RAD','rad_membcs',kfl_funno_rad)
     allocate(kfl_funbo_rad(nboun),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FUNBO_RAD','rad_membcs',kfl_funbo_rad)
     allocate(kfl_funty_rad(10),   stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FUNTY_RAD','rad_membcs',kfl_funty_rad)
     allocate(funpa_rad(6,10),     stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'FUNPA_RAD',    'rad_membcs',funpa_rad)

  case(3)
     !
     ! Deallocate memory
     !
     call memchk(two,istat,mem_modul(1:2,modul),'KFL_FIXNO_RAD','rad_membcs',kfl_fixno_rad)
     deallocate(kfl_fixno_rad,stat=istat)
     if(istat/=0) call memerr(two,'KFL_FIXNO_RAD','nsi_membcs',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'KFL_FIXBO_RAD','rad_membcs',kfl_fixbo_rad)
     deallocate(kfl_fixbo_rad,stat=istat)
     if(istat/=0) call memerr(two,'KFL_FIXBO_RAD','nsi_membcs',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'BVESS_RAD','rad_membcs',bvess_rad)
     deallocate(bvess_rad,stat=istat)
     if(istat/=0) call memerr(two,'BVESS_RAD','nsi_membcs',0_ip)
     if(nboun/=0) then
        call memchk(two,istat,mem_modul(1:2,modul),'BVNAT_RAD','rad_membcs',bvnat_rad)
        deallocate(bvnat_rad,stat=istat)
        if(istat/=0) call memerr(two,'BVNAT_RAD','rad_membcs',0_ip)
     end if
     if(kfl_conbc_rad==0)  then
        call memchk(two,istat,mem_modul(1:2,modul),'KFL_FUNNO_RAD','rad_membcs',kfl_funno_rad)
        deallocate(kfl_funno_rad,stat=istat)
        if(istat/=0) call memerr(two,'KFL_FUNNO_RAD','nsi_membcs',0_ip)
        call memchk(two,istat,mem_modul(1:2,modul),'KFL_FUNBO_RAD','rad_membcs',kfl_funbo_rad)
        deallocate(kfl_funbo_rad,stat=istat)
        if(istat/=0) call memerr(two,'KFL_FUNBO_RAD','nsi_membcs',0_ip)
        call memchk(two,istat,mem_modul(1:2,modul),'KFL_FUNTY_RAD','rad_membcs',kfl_funty_rad)
        deallocate(kfl_funty_rad,stat=istat)
        if(istat/=0) call memerr(two,'KFL_FUNTY_RAD','nsi_membcs',0_ip)
        call memchk(two,istat,mem_modul(1:2,modul),'FUNPA_RAD','rad_membcs',funpa_rad)
        deallocate(funpa_rad,stat=istat)
        if(istat/=0) call memerr(two,'FUNPA_RAD','nsi_membcs',0_ip)
     endif
        
  end select

end subroutine rad_membcs
