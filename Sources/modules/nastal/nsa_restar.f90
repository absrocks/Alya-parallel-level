subroutine nsa_restar(itask)
!------------------------------------------------------------------------
!****f* Nastal/nsa_restar
! NAME 
!    nsa_restar
! DESCRIPTION
!    This routine reads or writes the restart file
! USES
! USED BY
!    nsa_turnon
!***
!------------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_nastal
  use      mod_postpr
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iwopo,icomp,kfl_gores
  !
  ! Check if restrt file should be read or written
  !
  call respre(itask,kfl_gores)
  if( kfl_gores == 0 ) return
  !
  ! variables
  !
  icomp = min(TIME_N,ncomp_nsa)

  iwopo = 2

  gevec => umome(:,:,icomp)
  call postpr(gevec,postp(1) % wopos(1:2,iwopo),ittim,cutim)
  iwopo =  21
  gesca => densi(:,icomp)
  call postpr(gesca,postp(1) % wopos(1:2,iwopo),ittim,cutim)
  iwopo =  23
  gesca => energ(:,icomp)
  call postpr(gesca,postp(1) % wopos(1:2,iwopo),ittim,cutim)
  iwopo =  22  
  !
  ! Finish
  !
  call respre(3_ip,kfl_gores)

end subroutine nsa_restar
