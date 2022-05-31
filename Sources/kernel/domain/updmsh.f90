subroutine updmsh(itask)
!-----------------------------------------------------------------------
!****f* domain/updmsh
! NAME
!    updmsh
! DESCRIPTION
!    This routine performs update of mesh data. Depending on the 
!    refinement strategy, some data is saved (to be used in the
!    generation of the new mesh and projections) and some arrays
!    are deallocated (depending if they change or not)
! USES
! USED BY
!    Newmsh
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use mod_memchk
  implicit none
  integer(ip) :: itask,inode,inodb
  integer(ip) :: ielem,ielty,iboun
  integer(4)  :: istat 


end subroutine updmsh
