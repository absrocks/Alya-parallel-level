subroutine nsa_compline(x1,x2,y1,y2,m,q)
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_compline
  ! NAME 
  !    nsa_inimet
  ! DESCRIPTION
  !    This routine computes m and q for a segment 
  !    that passes through two points
  ! USED BY
  !    nsa_initial_conditiona
  !***
  !-----------------------------------------------------------------------
  use      def_master

  implicit none

  real(rp)    :: q,m,x1,x2,y1,y2
  
  !Compute m and q
  m = (y1 - y2)/(x1 - x2)
  q = y2 - m*x2

end subroutine nsa_compline
