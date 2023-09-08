subroutine rad_boumat(&
     pnode,pnodb,lboel,xmmat,xmrhs,gbsha,&
     gbsur,elmat,elrhs)
  !------------------------------------------------------------------------
  !****f* Radiat/rad_boumat
  ! NAME 
  !    rad_boumat
  ! DESCRIPTION
  !    Assemble boundary contribution
  ! USES
  ! USED BY
  !    rad_bouope
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  implicit none
  integer(ip), intent(in)    :: pnode,pnodb
  integer(ip), intent(in)    :: lboel(pnodb)
  real(rp),    intent(in)    :: xmmat,xmrhs,gbsha(pnodb),gbsur
  real(rp),    intent(inout) :: elmat(pnode,pnode),elrhs(pnode)
  integer(ip)                :: inodb,jnodb,inode,jnode
  real(rp)                   :: xmuit

  do inodb=1,pnodb  
     inode=lboel(inodb)
     elrhs(inode)=elrhs(inode)+gbsha(inodb)*xmrhs*gbsur
     xmuit=xmmat*gbsha(inodb)
     do jnodb=1,pnodb
        jnode=lboel(jnodb)
        elmat(jnode,inode)=elmat(jnode,inode)&
             +xmuit*gbsha(jnodb)*gbsur
     end do
  end do

end subroutine rad_boumat
