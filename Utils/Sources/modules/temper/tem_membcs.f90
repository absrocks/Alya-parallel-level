subroutine tem_membcs(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_membcs
  ! NAME
  !    tem_membcs
  ! DESCRIPTION
  !    Allocate memory for the physical problem
  ! OUTPUT 
  ! USES
  ! USED BY
  !    tem_reaphy
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_temper
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin
  integer(4)              :: istat

  select case(itask)

  case(1)

     if( INOTMASTER )then
        !
        ! Fixity and boundary values
        !
        allocate(kfl_fixno_tem(1,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXNO_TEM','tem_membcs',kfl_fixno_tem)
        do ipoin=1,npoin
           kfl_fixno_tem(1,ipoin)=-1
        end do
        allocate(kfl_fixbo_tem(nboun),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXBO_TEM','tem_membcs',kfl_fixbo_tem)
        if(kfl_conbc_tem==0) then
           allocate(bvess_tem(1,npoin,2),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_TEM','tem_membcs',bvess_tem)
           if(nboun/=0) then
              allocate(bvnat_tem(npnat_tem,nboun,2),stat=istat)
              call memchk(zero,istat,mem_modul(1:2,modul),'BVNAT_TEM','tem_membcs',bvnat_tem)
           end if
        else
           allocate(bvess_tem(1,npoin,1),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_TEM','tem_membcs',bvess_tem)
           if(nboun/=0) then
              allocate(bvnat_tem(npnat_tem,nboun,1),stat=istat)
              call memchk(zero,istat,mem_modul(1:2,modul),'BVNAT_TEM','tem_membcs',bvnat_tem)
           end if
        end if

     else
        ! 
        ! Allocate minimum size array if I am the master
        !
        allocate(kfl_fixno_tem(1,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXNO_TEM','tem_membcs',kfl_fixno_tem)
        allocate(kfl_fixbo_tem(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXBO_TEM','tem_membcs',kfl_fixbo_tem)
        if(kfl_conbc_tem==0) then
           allocate(bvess_tem(1,1,1),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_TEM','tem_membcs',bvess_tem)
           if(nboun/=0) then
              allocate(bvnat_tem(1,1,1),stat=istat)
              call memchk(zero,istat,mem_modul(1:2,modul),'BVNAT_TEM','tem_membcs',bvnat_tem)
           end if
        else
           allocate(bvess_tem(1,1,1),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_TEM','tem_membcs',bvess_tem)
           if(nboun/=0) then
              allocate(bvnat_tem(1,1,1),stat=istat)
              call memchk(zero,istat,mem_modul(1:2,modul),'BVNAT_TEM','tem_membcs',bvnat_tem)
           end if
        end if

     end if
        
  case(2)

     !
     ! Non-constant b.c.'s : Functions
     !
     allocate(kfl_funno_tem(npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FUNNO_TEM','tem_membcs',kfl_funno_tem)
     allocate(kfl_funbo_tem(nboun),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FUNBO_TEM','tem_membcs',kfl_funbo_tem)

  case(3)

     !
     ! Deallocate memory
     !
     call memchk(two,istat,mem_modul(1:2,modul),'KFL_FIXNO_TEM','tem_membcs',kfl_fixno_tem)
     deallocate(kfl_fixno_tem,stat=istat)
     if(istat/=0) call memerr(two,'KFL_FIXNO_TEM','nsi_membcs',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'KFL_FIXBO_TEM','tem_membcs',kfl_fixbo_tem)
     deallocate(kfl_fixbo_tem,stat=istat)
     if(istat/=0) call memerr(two,'KFL_FIXBO_TEM','nsi_membcs',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'BVESS_TEM','tem_membcs',bvess_tem)
     deallocate(bvess_tem,stat=istat)
     if(istat/=0) call memerr(two,'BVESS_TEM','nsi_membcs',0_ip)
     if(nboun/=0) then
        call memchk(two,istat,mem_modul(1:2,modul),'BVNAT_TEM','tem_membcs',bvnat_tem)
        deallocate(bvnat_tem,stat=istat)
        if(istat/=0) call memerr(two,'BVNAT_TEM','tem_membcs',0_ip)
     end if
        
  case(4)

     allocate(kfl_funty_tem(10),   stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FUNTY_TEM','tem_membcs',kfl_funty_tem)
     allocate(funpa_tem(6,10),     stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'FUNPA_TEM',    'tem_membcs',funpa_tem)

  end select

end subroutine tem_membcs
