subroutine got_elmma3(&
     pgaus,pnode,pevat,gpvol,resic,tescb,gprhs,elmat,elrhs)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmma3
  ! NAME 
  !    got_elmma3
  ! DESCRIPTION
  !    Assemble the elemental matrix from Gauss point contributions
  ! USES
  ! USED BY
  !    got_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)  :: pgaus,pnode,pevat
  real(rp),    intent(in)  :: gpvol(pgaus),gprhs(pgaus)
  real(rp),    intent(in)  :: resic(pnode,pgaus),tescb(pnode,pgaus)
  real(rp),    intent(out) :: elmat(pevat,pevat),elrhs(pevat)
  integer(ip)              :: igaus,inode,jnode
  real(rp)                 :: fact1
  !
  ! Initialization
  !
  do inode=1,pnode
     elrhs(inode)=0.0_rp
     do jnode=1,pnode
        elmat(jnode,inode)=0.0_rp
     end do
  end do
  !
  ! LHS and RHS
  !
  do igaus=1,pgaus 
     do inode=1,pnode
        fact1=gpvol(igaus)*tescb(inode,igaus)
        do jnode=1,pnode
           elmat(inode,jnode)=elmat(inode,jnode)&
                +fact1*resic(jnode,igaus)
        end do
        elrhs(inode)=elrhs(inode)+fact1*gprhs(igaus)       
     end do
  end do

end subroutine got_elmma3
