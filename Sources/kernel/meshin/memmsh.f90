subroutine memmsh(itask)
  !-----------------------------------------------------------------------
  !****f* meshin/memmsh
  ! NAME
  !    memmsh
  ! DESCRIPTION
  !    Allocate/Deallocate meshin arrays
  ! OUTPUT
  ! USED BY
  !    reamsh
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_domain, only    :  ndime
  use def_parame
  use def_master
  use def_inpout
  use def_meshin
  use mod_memchk
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case(itask)

  case(1_ip)
     !
     ! Sources
     !
     allocate(rsgeo(ndime,ndime,nsour),stat=istat)
     call memchk(zero,istat,memor_msh,'RSGEO','memmsh',rsgeo)              
     allocate(rsour(3,nsour),stat=istat)
     call memchk(zero,istat,memor_msh,'RSOUR','memmsh',rsour)

  case(2_ip)
     !
     ! Allocate memory for renum
     !

  end select

end subroutine memmsh
