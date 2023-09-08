subroutine par_bridge(itask)
  !-----------------------------------------------------------------------
  !****f* Parall/par_bridge
  ! NAME
  !    par_bridge
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    Reapro
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_parall
  use def_inpout
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case(itask)

  case(1_ip)
     !
     ! Allocate memory and slaves receive
     !
     if( parii == 1 ) then
        allocate(parin(npari),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,5),'parin','par_bridge',parin)
        allocate(parre(nparr),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,5),'parre','par_bridge',parre)
        if( ISLAVE .or. kfl_ptask==2 ) call par_broadc()
     end if

  case(2_ip)
     !
     ! Master broadcast and deallocate memory
     !
     if( IMASTER .and. kfl_ptask/=2 ) call par_broadc()
     call memchk(two,istat,mem_servi(1:2,5),'parin','par_bridge',parin)
     deallocate(parin,stat=istat)
     if( istat /= 0 ) call memerr(two,'parin','par_bridge',0_ip)
     call memchk(two,istat,mem_servi(1:2,5),'parre','par_bridge',parre)
     deallocate(parre,stat=istat)
     if( istat /= 0 ) call memerr(two,'parre','par_bridge',0_ip)     

  end select

end subroutine par_bridge
