subroutine elsest_elbbox(&
     mnode,ndime,npoin,nelem,nnode,lnods,ltype,coord,xmima)
  !
  ! Construct the bin structure
  !
  use def_elsest
  implicit none
  integer(ip), intent(in)  :: mnode,ndime,npoin,nelem
  integer(ip), intent(in)  :: nnode(*)
  integer(ip), intent(in)  :: lnods(mnode,*),ltype(*)
  real(rp),    intent(in)  :: coord(ndime,npoin)
  real(rp),    intent(out) :: xmima(3,2,nelem)
  integer(ip)              :: ielem,inode,ipoin
  real(rp)                 :: zeror

  zeror = epsilon(1.0_rp)

  if( ndime == 1 ) then

     do ielem = 1,nelem
        xmima(1,1,ielem) = coord(1,lnods(1,ielem))
        xmima(1,2,ielem) = coord(1,lnods(1,ielem))
        do inode = 2,nnode(abs(ltype(ielem)))
           ipoin            = lnods(inode,ielem)
           xmima(1,1,ielem) = min( xmima(1,1,ielem) , coord(1,ipoin) )
           xmima(1,2,ielem) = max( xmima(1,2,ielem) , coord(1,ipoin) )
        end do
        xmima(1,1,ielem) = xmima(1,1,ielem) - zeror
        xmima(1,2,ielem) = xmima(1,2,ielem) + zeror
     end do

  else if( ndime == 2 ) then

     do ielem = 1,nelem
        xmima(1,1,ielem) = coord(1,lnods(1,ielem)) 
        xmima(2,1,ielem) = coord(2,lnods(1,ielem))
        xmima(1,2,ielem) = coord(1,lnods(1,ielem))
        xmima(2,2,ielem) = coord(2,lnods(1,ielem))
        do inode = 2,nnode(abs(ltype(ielem)))
           ipoin            = lnods(inode,ielem)
           xmima(1,1,ielem) = min( xmima(1,1,ielem) , coord(1,ipoin) )
           xmima(2,1,ielem) = min( xmima(2,1,ielem) , coord(2,ipoin) )
           xmima(1,2,ielem) = max( xmima(1,2,ielem) , coord(1,ipoin) )
           xmima(2,2,ielem) = max( xmima(2,2,ielem) , coord(2,ipoin) )
        end do
        xmima(1,1,ielem) = xmima(1,1,ielem) - zeror
        xmima(1,2,ielem) = xmima(1,2,ielem) + zeror
        xmima(2,1,ielem) = xmima(2,1,ielem) - zeror
        xmima(2,2,ielem) = xmima(2,2,ielem) + zeror
     end do

  else

     do ielem = 1,nelem
        xmima(1,1,ielem) = coord(1,lnods(1,ielem)) 
        xmima(2,1,ielem) = coord(2,lnods(1,ielem))  
        xmima(3,1,ielem) = coord(3,lnods(1,ielem))  
        xmima(1,2,ielem) = coord(1,lnods(1,ielem)) 
        xmima(2,2,ielem) = coord(2,lnods(1,ielem)) 
        xmima(3,2,ielem) = coord(3,lnods(1,ielem)) 
        do inode = 2,nnode(abs(ltype(ielem)))
           ipoin            = lnods(inode,ielem)
           xmima(1,1,ielem) = min( xmima(1,1,ielem) , coord(1,ipoin) )
           xmima(2,1,ielem) = min( xmima(2,1,ielem) , coord(2,ipoin) )
           xmima(3,1,ielem) = min( xmima(3,1,ielem) , coord(3,ipoin) )
           xmima(1,2,ielem) = max( xmima(1,2,ielem) , coord(1,ipoin) )
           xmima(2,2,ielem) = max( xmima(2,2,ielem) , coord(2,ipoin) )
           xmima(3,2,ielem) = max( xmima(3,2,ielem) , coord(3,ipoin) )
        end do
        xmima(1,1,ielem) = xmima(1,1,ielem) - zeror
        xmima(1,2,ielem) = xmima(1,2,ielem) + zeror
        xmima(2,1,ielem) = xmima(2,1,ielem) - zeror
        xmima(2,2,ielem) = xmima(2,2,ielem) + zeror
        xmima(3,1,ielem) = xmima(3,1,ielem) - zeror
        xmima(3,2,ielem) = xmima(3,2,ielem) + zeror
     end do
    
  end if
end subroutine elsest_elbbox
