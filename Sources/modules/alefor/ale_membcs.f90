subroutine ale_membcs(itask)
  !-----------------------------------------------------------------------
  !****f* Domain/ale_membcs
  ! NAME
  !    ale_membcs
  ! DESCRIPTION
  !    Allocate/Deallocate the geometry arrays 
  !    ITASK=1 ... Allocate memory
  !    ITASK=2 ... Deallocate memory
  ! OUTPUT
  ! USED BY
  !    reageo
  !    sengeo
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  use mod_memchk
  use mod_memory
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ifunc
  integer(4)              :: istat

  select case(itask)

  case(  1_ip )
     !
     ! Fixity and boundary values
     ! 
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_ALE','ale_membcs',kfl_fixno_ale,ndime,max(1_ip,npoin))
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXBO_ALE','ale_membcs',kfl_fixbo_ale,max(1_ip,nboun))
     call memory_alloca(mem_modul(1:2,modul),'BVESS_ALE'    ,'ale_membcs',bvess_ale,ndime,max(1_ip,npoin))

  case( 5_ip)
     !
     ! Support geometry
     !
     call memory_alloca(mem_modul(1:2,modul),'COORD_AD','ale_membcs' , coord_ad , ndime    , npoin_ad )
     call memory_alloca(mem_modul(1:2,modul),'LNODB_AD','ale_membcs' , lnodb_ad , mnodb_ad , nboun_ad )
     call memory_alloca(mem_modul(1:2,modul),'LTYPB_AD','ale_membcs' , ltypb_ad , nboun_ad )

  case(-5_ip)
     !
     ! Support geometry
     !
     call memory_deallo(mem_modul(1:2,modul),'COORD_AD','ale_membcs' , coord_ad )
     call memory_deallo(mem_modul(1:2,modul),'LNODB_AD','ale_membcs' , lnodb_ad )
     call memory_deallo(mem_modul(1:2,modul),'LTYPB_AD','ale_membcs' , ltypb_ad )

  case( 38_ip)
     !
     ! KFL_FIXBO_ALE: Allocate kfl_fixbo_ale
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXBO_ALE','ale_membcs',kfl_fixbo_ale,max(1_ip,nboun))

  case(-38_ip)
     !
     ! KFL_FIXBO_ALE: Dellocate kfl_fixbo_ale
     !
     call memory_deallo(mem_modul(1:2,modul),'KFL_FIXBO_ALE','ale_membcs',kfl_fixbo_ale)

  case( 44_ip)
     !
     ! BVESS_ALE: Allocate bvess_ale
     !
     allocate(bvess_ale(1,npoin),stat=istat)
     call memchk(zero,istat,memor_dom,'BVESS_ALE','ale_membcs',bvess_ale)

  case(-44_ip)
     !
     ! BVESS_ALE: Dellocate bvess_ale
     !
     deallocate(bvess_ale,stat=istat)
     if(istat/=0) call memerr(two,'BVESS_ALE','ale_membcs',0_ip)
     call memchk(two,istat,memor_dom,'BVESS_ALE','ale_membcs',bvess_ale)

  case( 45_ip)
     !
     ! FUNNO_ALE: Allocate funno_ale
     !
     allocate(kfl_funno_ale(npoin),stat=istat)
     call memchk(zero,istat,memor_dom,'KFL_FUNNO_ALE','ale_membcs',kfl_funno_ale)

  case(-45_ip)
     !
     ! FUNNO_ALE: Dellocate funno_ale
     !
     deallocate(kfl_funno_ale,stat=istat)
     if(istat/=0) call memerr(two,'KFL_FUNNO_ALE','ale_membcs',0_ip)
     call memchk(two,istat,memor_dom,'KFL_FUNNO_ALE','ale_membcs',kfl_funno_ale)

  case(46_ip)
     !
     ! FUNPA_ALE
     !
     ifunc = igene
     if( kfl_funty_ale(ifunc,1) /= 0 ) then   
        allocate(funpa_ale(ifunc)%a(kfl_funty_ale(ifunc,2)),stat=istat)
        call memchk(zero,istat,memor_dom,'FUNPA_ALE','ale_membcs',funpa_ale(ifunc)%a) 
     end if

  case(-46_ip)
     !
     ! FUNPA_ALE 
     !
     ifunc = igene
     if( kfl_funty_ale(ifunc,1) /= 0 ) then 
        call memchk(two,istat,memor_dom,'funpa_ale(ifunc)%a','ale_membcs',funpa_ale(ifunc)%a)
        deallocate(funpa_ale(ifunc)%a,stat=istat)
        if( istat /= 0 ) call memerr(two,'funpa_ale(ifunc)%a','ale_membcs',0_ip)
     end if

  case(47_ip)
     !
     ! Time/space functions
     !
     allocate(kfl_funty_ale(mfunc_ale,2), stat=istat)
     call memchk(zero,istat,memor_dom,'KFL_FUNTY_ALE','ale_membcs',kfl_funty_ale)
     allocate(funpa_ale(mfunc_ale),stat=istat)
     call memchk(zero,istat,memor_dom,'FUNPA_ALE',    'ale_membcs',funpa_ale)

  case(-47_ip)
     !
     ! Time/space functions
     !
     call memchk(two,istat,memor_dom,'KFL_FUNTY_ALE','ale_membcs',kfl_funty_ale)
     deallocate(kfl_funty_ale,stat=istat)
     if( istat /= 0 ) call memerr(two,'KFL_FUNTY_ALE','ale_membcs',0_ip)
     call memchk(two,istat,memor_dom,'FUNPA_ALE','ale_membcs',funpa_ale)
     deallocate(funpa_ale,stat=istat)
     if( istat /= 0 ) call memerr(two,'FUNPA_ALE','ale_membcs',0_ip)

  case( 48_ip)
     !
     ! FUNBO_ALE: Allocate funbo_ale
     !
     allocate(kfl_funbo_ale(nboun),stat=istat)
     call memchk(zero,istat,memor_dom,'KFL_FUNBO_ALE','ale_membcs',kfl_funbo_ale)

  case(-48_ip)
     !
     ! FUNBO_ALE: Dellocate funbo_ale
     !
     deallocate(kfl_funbo_ale,stat=istat)
     if(istat/=0) call memerr(two,'KFL_FUNBO_ALE','ale_membcs',0_ip)
     call memchk(two,istat,memor_dom,'KFL_FUNBO_ALE','ale_membcs',kfl_funbo_ale)

  end select

end subroutine ale_membcs
