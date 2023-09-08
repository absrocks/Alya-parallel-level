subroutine elmdir(&
     pnode,lnods,kfl_fixno,bvess,elmat,elrhs)
  !------------------------------------------------------------------------
  !****f* mathru/elmdir
  ! NAME 
  !    elmdir
  ! DESCRIPTION
  !    This routine prescribes the boundary conditions for the 
  !    temperature equations. 
  ! USES
  ! USED BY
  !    elmadr
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  implicit none
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: lnods(pnode)
  integer(ip), intent(in)    :: kfl_fixno(*)
  real(rp),    intent(in)    :: bvess(*)
  real(rp),    intent(inout) :: elmat(pnode,pnode),elrhs(pnode)
  integer(ip)                :: inode,ipoin,jnode
  real(rp)                   :: adiag
  
  do inode=1,pnode
     ipoin=lnods(inode)
     if(kfl_fixno(ipoin)>=1) then
        adiag=elmat(inode,inode)
        do jnode=1,pnode
           elmat(inode,jnode)=0.0_rp
           elrhs(jnode)=elrhs(jnode)&
                -elmat(jnode,inode)*bvess(ipoin)
           elmat(jnode,inode)=0.0_rp
        end do
        elmat(inode,inode)=adiag
        elrhs(inode)=adiag*bvess(ipoin)
     end if
  end do

end subroutine elmdir
