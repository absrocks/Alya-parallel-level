subroutine setlbe
  !-----------------------------------------------------------------------
  !****f* Domain/setlbe
  ! NAME
  !    setlbe
  ! DESCRIPTION
  !    This routine constucts boundary arrays
  ! OUTPUT
  !    LBOEL(MNODB,NBOUN) ... Boundary/Element connectivity
  !    LELBO(NBOUN) ......... Boundary/ element 
  !    LTYPB(NBOUN) ......... Type of boundary
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  implicit none
  integer(ip)          :: iboun,ielem,inodb,ipoin,jpoin,inode,idime
  integer(ip)          :: knodb,knode

  real(rp)             :: disti

  if( INOTMASTER ) then     
     ! 
     ! LBOEL if the elements connected to the boundaries are unknown
     ! This does not work with internal boundaries as lbouel could
     ! pick up the element on the wrong side
     !
     !!!!!!!call lbouel()
     !
     ! Calculate only the boundary/element nodal connectivity
     !
     call memgeo(5_ip)
     do iboun = 1,nboun
        knodb = nnode(ltypb(iboun))
        ielem = lelbo(iboun)
        knode = nnode(ltype(ielem))
        do inodb = 1,knodb
           ipoin = lnodb(inodb,iboun)
           nodes3: do inode = 1,knode
              jpoin = lnods(inode,ielem)
              if( ipoin == jpoin ) then
                 lboel(inodb,iboun) = inode
                 exit nodes3
              end if
           end do nodes3
        end do
     end do
     !
     ! Order boundaries
     !
     call order_boundaries()
     !
     ! Check element/boundary local numbering
     !
     call mescek(7_ip)

  end if

100 format(20(2x,i12))
200 format(i12,2x,20(2x,f15.8))

end subroutine setlbe
