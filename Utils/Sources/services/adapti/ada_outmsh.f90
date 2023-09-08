subroutine ada_outmsh
!-----------------------------------------------------------------------
!****f* adapti/ada_outmsh
! NAME 
!    ada_outmsh
! DESCRIPTION
!    This routine outputs the modified dom.dat file
!    
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

  if (nimmo_ada > 0) then
     !
     ! Immersed objects are present
     !
     call ada_outmgi
  end if

  call ada_outgid       ! gid post files, just to check everything went right
  call ada_outaly       ! alya dom.geo and fix files

end subroutine ada_outmsh
