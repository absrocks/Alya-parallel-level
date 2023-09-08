subroutine rad_elmdir(&
     pnode,lnods,elmat,elrhs)
  !------------------------------------------------------------------------
  !****f* Radiat/rad_elmdir
  ! NAME 
  !    rad_elmdir
  ! DESCRIPTION
  !    This routine prescribes Dirichlet boundary conditions for the 
  !    radiation heat transfer equations (only used for exact solutions)
  ! USES
  ! USED BY
  !    rad_elmope
  !    rad_bouope
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_radiat, only       :  bvess_rad,kfl_fixno_rad
  implicit none
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(inout) :: elmat(pnode,pnode),elrhs(pnode)
  integer(ip)                :: inode,ipoin,jnode
  real(rp)                   :: adiag,xvalu
  
!!$  do inode = 1,pnode
!!$     ipoin = lnods(inode)
!!$     if(  kfl_fixno_rad(1,ipoin) == 1 .or.&
!!$          kfl_fixno_rad(1,ipoin) == 4 .or.&
!!$          kfl_fixno_rad(1,ipoin) == 5 ) then
!!$        adiag = elmat(inode,inode)
!!$        xvalu = bvess_rad(ipoin,1)
!!$        do jnode = 1,pnode
!!$           elmat(inode,jnode) = 0.0_rp
!!$           elrhs(jnode)       = elrhs(jnode) - elmat(jnode,inode) * xvalu
!!$           elmat(jnode,inode) = 0.0_rp
!!$        end do
!!$        elmat(inode,inode) = adiag
!!$        elrhs(inode)       = adiag * xvalu
!!$     end if
!!$  end do

end subroutine rad_elmdir
