subroutine elmmas(pnode,pgaus,gpvol,gpsha,elmat)
  !-----------------------------------------------------------------------
  !****f* domain/elmmas
  ! NAME
  !    elmmas
  ! DESCRIPTION
  !    This routines calculates the element mass matrix (symmetric) 
  ! OUTPUT
  !    ELMAT : Consistent mass matrix
  ! USED BY
  !    Domain
  !*** 
  !-----------------------------------------------------------------------
  use def_kintyp, only      : rp,ip
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus
  real(rp),    intent(in)  :: gpvol(pgaus),gpsha(pnode,pgaus)
  real(rp),    intent(out) :: elmat(pnode,pnode)
  integer(ip)              :: igaus,inode,jnode
  real(rp)                 :: fact1

  do inode=1,pnode
     do jnode=1,pnode
        elmat(inode,jnode)=0.0_rp
     end do
  end do

  do igaus=1,pgaus
     do inode=1,pnode
        fact1=gpvol(igaus)*gpsha(inode,igaus)
        do jnode=1,pnode
           elmat(inode,jnode)=elmat(inode,jnode)&
                +fact1*gpsha(jnode,igaus)
        end do     
     end do
  end do

end subroutine elmmas
