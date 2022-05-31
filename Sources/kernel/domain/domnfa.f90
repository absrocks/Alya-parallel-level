subroutine domnfa(ndime,nnode,nface)
!-----------------------------------------------------------------------
!****f* Domain/domnfa
! NAME
!    domnfa
! DESCRIPTION
!    Compute number of faces of an element
! OUTPUT
!    NFACE : # of faces of an element
! USED BY
!    Domain
!***
!-----------------------------------------------------------------------
  use      def_kintyp
  implicit none 
  integer(ip), intent(in)  ::  ndime,nnode
  integer(ip), intent(out) ::  nface
      
  if(ndime==2) then 
     !
     ! 2D elements
     !
     if (nnode<=4) then 
        nface=nnode
     else
        if (nnode==6.or.nnode==7) then
           nface=3
        else
           nface=4
        end if
     end if
  else
     !
     ! 3D elements
     !
     if (nnode==4) then
        nface=4
     else if (nnode== 5) then
        nface=5
     else if (nnode==10) then
        nface=4
     else if (nnode==18) then
        nface=5
     else if (nnode==8) then
        nface=6
     else if (nnode==20) then
        nface=6
     else if (nnode==27) then
        nface=6
     else if (nnode==6) then
        nface=5
     end if
  end if

end subroutine domnfa

