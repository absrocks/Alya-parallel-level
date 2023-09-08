subroutine got_memphy(itask)
!-----------------------------------------------------------------------
!****f* Gotita/got_memphy
! NAME 
!    got_memphy
! DESCRIPTION
!    This routine allocates memory for physical arrays
! USES
!    ecoute
!    memchk
!    runend
! USED BY
!    got_reaphy
!***
!-----------------------------------------------------------------------
  use def_parame 
  use def_inpout
  use def_master
  use def_gotita
  use def_domain
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case(itask)

  case(1)
     allocate(veloc_got(ndime,npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VELOC_GOT','got_memphy',veloc_got)

  case(2)
     allocate(vdrop(ndime,npoin,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VDROP','got_memphy',vdrop)

  end select

end subroutine got_memphy
