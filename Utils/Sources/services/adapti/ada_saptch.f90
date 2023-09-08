subroutine ada_saptch
!-----------------------------------------------------------------------
!****f* adapti/ada_saptch
! NAME 
!    ada_saptch
! DESCRIPTION
!    This routine fills the hole left by ada_sahole with a patch mesh
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
  integer(ip)   :: &
       ielem,ipoin,ielsu,ielad,inode,jnode,ifono,kfono
  integer(ip)   :: &
       lnofa(mnode)



  !
  ! Make the announcement of success
  !
  call ada_livinf('mesh adapted... leaving!')
  

end subroutine ada_saptch
