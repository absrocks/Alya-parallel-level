subroutine reasla
!-----------------------------------------------------------------------
!****f* Domain/reasla
! NAME
!    reasla
! DESCRIPTION
!    LMASL(IPOIN) = IMASL  >  0: IPOIN is the master of IMASL
!                          <  0: IPOIN is the slave of -IMASL
!                          =  0: IPOIN is neither a master nor a slave
!    This routine also calls to the routines that creare R_SOL, C_SOL
!    and PERMX.
! USED BY
!    cresla
!***
!-----------------------------------------------------------------------
  use      def_domain
  use      def_master
  use      def_inpout
  implicit none
!  integer(ip) :: islav,imast
!  real(rp)    :: vslav(9),vmast(9)
!  logical(lg) :: slerr
      
end subroutine reasla

