subroutine exm_outrep
!-----------------------------------------------------------------------
!****f* Exmedi/exm_outrep
! NAME 
!    exm_outrep
! DESCRIPTION
!    This routine write reports to the report files.
! USES
!
! USED BY
!    exm_output (itask=0,1)
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master

  use      def_exmedi

  implicit none
  real(rp)    :: xvalu(200)
  integer(ip) :: ivalu(200),nrepi,nreva,ireva,ireco,kreco, & 
       nreco,irepi,ipoin,kdrep

end subroutine exm_outrep
