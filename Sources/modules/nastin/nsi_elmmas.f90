subroutine nsi_elmmas(&
     pnode,pgaus,gpsha,gpvol,gpden,elmas)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmma4
  ! NAME 
  !    nsi_elmma4
  ! DESCRIPTION
  !    Compute element matrix and RHS
  ! USES
  ! USED BY
  !    nsi_elmop3
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only       :  ip,rp
  implicit none
  integer(ip), intent(in)    :: pnode,pgaus
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(out)   :: elmas(pnode)
  integer(ip)                :: inode,igaus
  real(rp)                   :: xfact

     elmas = 0.0_rp
     do igaus = 1,pgaus
        xfact = gpvol(igaus) * gpden(igaus)
        do inode = 1,pnode
           elmas(inode) = elmas(inode) + gpsha(inode,igaus) * xfact
        end do
     end do

end subroutine nsi_elmmas
