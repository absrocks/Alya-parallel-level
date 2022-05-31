subroutine elmlum(pnode,pgaus,pelty,elcod,gprea,gpcar,elmat)
  !-----------------------------------------------------------------------
  !****f* domain/elmlum
  ! NAME
  !    elmlum
  ! DESCRIPTION
  !    This routines calculates the reaction term in element matrix
  !    using lumped mass. Reaction is first extrapolated from nodes to 
  !    nodes.
  ! OUTPUT
  !    ELMAT: Consistent mass matrix
  ! USED BY
  !    Domain
  !*** 
  !-----------------------------------------------------------------------
  use def_kintyp, only      : rp,ip
  use def_domain, only      : ndime,elmar
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus,pelty
  real(rp),    intent(in)  :: elcod(ndime,pnode)
  real(rp),    intent(in)  :: gprea(pgaus)
  real(rp),    intent(out) :: gpcar(ndime,pnode) 
  real(rp),    intent(out) :: elmat(pnode,pnode)
  integer(ip)              :: igaus,inode
  real(rp)                 :: xjaci(9),xjacm(9),react,gpvol

  do inode = 1,pnode
     react = 0.0_rp
     do igaus = 1,pgaus
        react = react + gprea(igaus) * elmar(pelty)%shaga(igaus,inode)
     end do
     call elmder(&
          pnode,ndime,elmar(pelty)%deric(1,1,inode),&
          elcod,gpcar,gpvol,xjacm,xjaci)
     gpvol = elmar(pelty)%weigc(inode)*gpvol
     elmat(inode,inode) = elmat(inode,inode) + react * gpvol
  end do

end subroutine elmlum
