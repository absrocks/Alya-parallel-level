 subroutine tem_elmext(pnode,elmat,elrhs)
  !-----------------------------------------------------------------------
  !****f* Nastin/tem_elmext
  ! NAME 
  !    tem_elmext
  ! DESCRIPTION
  !    MOdify element matrix when extension elements are used
  ! USES
  ! USED BY
  !    tem_elmext
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  implicit none
  integer(ip), intent(in)    :: pnode
  real(rp),    intent(out)   :: elmat(pnode,pnode)
  real(rp),    intent(out)   :: elrhs(pnode)
  integer(ip)                :: inode,jnode,idofn,jdofn
  !
  ! A, b
  !
  do inode = 2,pnode
     do jnode = 1,pnode
        elmat(inode,jnode) = 0.0_rp
     end do
     elrhs(inode) = 0.0_rp
  end do
  
end subroutine tem_elmext
