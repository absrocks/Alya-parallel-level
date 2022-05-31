subroutine assric(ndofn,pnode,lnods,elrhs,elmat,elunk,rhsid)
  !------------------------------------------------------------------------
  !****f* mathru/assric
  ! NAME 
  !    Assemble residual for matrix-free richardson solver
  ! DESCRIPTION
  !    Assembly of the RHS
  ! USES
  ! USED BY
  !    *_elmope
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only        :  ip,rp 
  use def_domain, only        :  npoin
  implicit none
  integer(ip),  intent(in)    :: ndofn,pnode
  integer(ip),  intent(in)    :: lnods(pnode)
  real(rp),     intent(inout) :: elrhs(*)
  real(rp),     intent(in)    :: elmat(pnode*ndofn,pnode*ndofn)
  real(rp),     intent(in)    :: elunk(*)
  real(rp),     intent(inout) :: rhsid(*)
  integer(ip)                 :: inode,jnode,ipoin,idofl,idofg,idofn,ndof2

  if(ndofn==1) then
     !
     ! 1 DOF
     !
     do inode=1,pnode
        ipoin=lnods(inode)
        do jnode=1,pnode
           elrhs(inode) = elrhs(inode) - elmat(inode,jnode) * elunk(jnode)
        end do
        rhsid(ipoin) = rhsid(ipoin) + elrhs(inode)
     end do
  else
     call runend('ASSRIC: NOT CODED')
  end if

end subroutine assric
