subroutine gtband(npoin,c_dom,r_dom,bandw,profi)
  !-----------------------------------------------------------------------
  !****f* Domain/gtband
  ! NAME
  !    gtband
  ! DESCRIPTION
  !    This subroutine computed the bandwth and profile of the mesh
  ! OUTPUT
  !    BANDW
  !    PROFI
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  implicit none
  integer(ip), intent(in)    :: npoin
  integer(ip), intent(in)    :: c_dom(*),r_dom(npoin+1)
  integer(ip), intent(inout) :: bandw
  real(rp),    intent(inout) :: profi
  integer(ip)                :: ipoin,izdomin,band,izdom,bandloc,ipmax,jpmax  

  bandw =  0_ip
  profi =  0.0_rp
  !naved =  0      ! Average number of edges
  !nmied =  1e6    ! Min number of edges
  !nmaed = -1e6    ! Max number of edges

  do ipoin = 1,npoin
     !
     ! Initialize local bandwth
     ! 
     bandloc = 0_ip
     !
     ! Loop on neighbors
     !
     do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
        izdomin = c_dom(izdom)
        if( ipoin /= izdomin ) then
           band = abs(izdomin-ipoin)
           !
           ! Test bandwth
           !
           if( band > bandw ) then
              bandw = band
              ipmax = ipoin
              jpmax = izdomin
           endif
           !
           ! Test profile
           !
           if( izdomin < ipoin ) then
              if( band > bandloc ) bandloc = band
           end if
        end if
     end do
     !
     ! Accumulate profile
     !
     profi = profi + real(bandloc)

  end do

end subroutine gtband
