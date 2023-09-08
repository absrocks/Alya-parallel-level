subroutine nsi_parall(itask)
!-----------------------------------------------------------------------
!****f* Nastin/nsi_parall
! NAME
!    nsi_parall
! DESCRIPTION
!    This routine is a bridge to Parall service  
! USED BY
!    Nastin
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_nastin
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipart,istar,istop

  if( IPARALL ) then

     select case (itask)
        
     case ( 1_ip )
        !
        ! Receive/Read data in nsi_reaphy, nsi_reanut and nsi_reaous
        !
        call nsi_sendat(one)

     end select
     
  end if

end subroutine nsi_parall
