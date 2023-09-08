subroutine ada_scagec(ipoin,iele1,iimmo,locsh)
!-----------------------------------------------------------------------
!****f* adapti/ada_scagec
! NAME 
!    ada_modmsh
! DESCRIPTION
!    This routine scatters geometrical connections
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain

  use      def_adapti
  implicit none

  integer(ip)    :: ipoin,iele1,inod1,iimmo
  real(rp)       :: locsh(mnode),xaux

  do inod1= 1,mnode
     gshac_ada(inod1,ipoin,iimmo)=locsh(inod1)
     gpoic_ada(inod1,ipoin,iimmo)=lnods(inod1,iele1)
  end do


end subroutine ada_scagec
