subroutine chm_elmdir(iclas,pnode,lnods,elmat,elrhs)
  !-----------------------------------------------------------------------
  !****f* Turbul/chm_elmdir
  ! NAME
  !   chm_elmdir
  ! DESCRIPTION
  ! This routine prescribes the boundary conditions for the 
  ! scalar equations. 
  ! USES
  ! USED BY
  !    chm_elmop1
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_chemic, only       :  kfl_fixno_chm,bvess_chm
  implicit none
  integer(ip), intent(in)    :: iclas,pnode
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(inout) :: elmat(pnode,pnode),elrhs(pnode)
  integer(ip)                :: inode,ipoin,jnode
  real(rp)                   :: adiag,xvalu

  do inode = 1,pnode
     ipoin = lnods(inode)
     if(  kfl_fixno_chm(iclas,ipoin)>=1) then
        xvalu = bvess_chm(iclas,ipoin)
        adiag = elmat(inode,inode)
        do jnode = 1,pnode
           elmat(inode,jnode)  = 0.0_rp
           elrhs(jnode)        = elrhs(jnode)-elmat(jnode,inode)*xvalu
           elmat(jnode,inode)  = 0.0_rp
        end do
        elmat(inode,inode) = adiag
        elrhs(inode)       = adiag*xvalu
     end if
  end do

end subroutine chm_elmdir
